################################################################################
##                Fit cumulative risk model to simulated data                 ##
################################################################################

library(sandwich)
library(survey)
library(kableExtra)
pal <- c("#F6AD4F", "#6CB9A9", "#49A5D6", "#A45336")

# Set your working directory
wd <- ''

# Load some useful functions used during the estimation approach
source(file = paste0(wd, 'functions_ouf_cr.R'))
# Source C++ code for calculating OU process precision matrix
Rcpp::sourceCpp(paste0(wd, "ou_precision_matrix.cpp"))

# Select set of true parameters under which data were simulated
setting <- 1 # OU process has higher correlation
# setting <- 2 # OU process has lower correlation

# Set seed
g <- 1
set.seed(412+g)

################################################################################
## Read in data
################################################################################

# dat_long contains info on the longitudinal outcomes
dat_long <- read.csv(paste0(wd, 'data/sim_dat_long_setting', setting, '.csv'))
# - user.id = ID for each individual
# - time = time of measurement occasion
# - active-scared = measured longitudinal outcomes

# Set sample size (number of independent individuals)
N <- length(unique(dat_long$user.id))

emotion_items <- c('active', 'calm', 'determined', 'enthusiastic', 'grateful',
                   'happy', 'proud', 'joyful', 'attentive',
                   'angry', 'ashamed', 'disgusted', 'guilty', 'irritable',
                   'lonely', 'nervous', 'sad', 'scared')

################################################################################
## Stage 1 estimation: fit OU factor model
################################################################################

# Parameters can either be initialized empirically or with user-defined values.
#  Specify the approach to be used here:
init_method <- 'empirical' # options are 'empirical' or 'user_def'
k <- 18 # number of longitudinal outcomes
p <- 2 # number of latent factors

source(paste0(wd, 'fit_ouf.R'))

stage1_ests <- res %>% filter(faboup_iter == max(res$faboup_iter))
# format parameter estimates from longitudinal submodel for predicting factors scores
theta_boup_vec <- stage1_ests %>%
  filter(variable %in% paste0(paste0('theta_ou_', 1:p), sort(rep(1:p, p))))
sigma_boup_vec <- stage1_ests %>%
  filter(variable %in% paste0(paste0('sigma_ou_', 1:p), 1:p))
lambda_vec <- stage1_ests %>%
  filter(variable %in% paste0('lambda_', 1:k))
sigma2_u_vec <- stage1_ests %>%
  filter(variable %in% paste0('sigma2_u_', 1:k))
sigma2_e_vec <- stage1_ests %>%
  filter(variable %in% paste0('sigma2_e_', 1:k))

# these are the final estimates from stage 1
# OU drift parameter
theta_boup <- matrix(theta_boup_vec$est, p, p)
# OU volatility parameter
sigma_boup <- diag(sigma_boup_vec$est)
# Loadings matrix
Lambda <- matrix(c(lambda_vec$est[1:9], rep(0, k),
                   lambda_vec$est[10:18]), ncol = p)
# Covariance matrix for random intercepts
Sigma_u <- diag(sigma2_u_vec$est)
# Covariance matrix for measurement error
Sigma_e <- diag(sigma2_e_vec$est)

################################################################################
## Predict factor scores at measurement occasions
################################################################################

# Calculate number of measurement occasions (for long. outcomes) per individual
measurement_n <- as.numeric(table(dat_long$user.id))
max_ni <- max(measurement_n)

# Reformat measurement times for each individual
measurement_times <- matrix(nrow = N, ncol = max_ni)
for (i in 1:N){
  ni <- measurement_n[i]
  measurement_times[i,1:ni] <- dat_long$time[dat_long$user.id == i]
}

# Calculate covariance matrix for latent factors (OU process) for each individual
Psi_list <- lapply(1:N, FUN = function(i){
  ni <- measurement_n[i]
  meas_times_i <- measurement_times[i, 1:ni]
  theta_t = t(theta_boup)
  sigma2_vec = matrix(sigma_boup %*% t(sigma_boup), ncol = 1)
  kron_sum_theta = kronsum(theta_boup, theta_boup)
  ni = length(c(meas_times_i))
  Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta, theta = theta_boup,
                                theta_t = theta_t, sigma2_vec = sigma2_vec,
                                ni = ni, times = c(meas_times_i))
  Psi_i <- solve(Omega_i)
  Psi_i <- (Psi_i + t(Psi_i))/2
  Psi_i # covariance matrix for OU process
})

# Calculate covariance matrix for observed longitudinal outcome for each
#  individual using OU factor model parameters
CovX_list <- lapply(1:N, FUN = function(i){
  meas_times_i <- dat_long$time[which(dat_long$user.id == i)]
  ni <- length(meas_times_i)
  Psi_i <- Psi_list[[i]] # covariance matrix for latent factors
  Lambda_mat <- kronecker(diag(ni), Lambda) %*% Psi_i %*% kronecker(diag(ni), t(Lambda))
  Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), Sigma_u)
  Sigma_e_mat <- kronecker(diag(ni), Sigma_e)
  
  # Cov(X_i) is the covariance for observed long. outcomes for individual i
  Lambda_mat + Sigma_u_mat + Sigma_e_mat # = Cov(X_i)
})

# Set up matrix of observed longitudinal outcomes
all_person_X <- matrix(NA, nrow = N, ncol = k*max_ni)
for (i in 1:N){
  obs_X <- dat_long %>%
    filter(user.id == i) %>%
    dplyr::select(all_of(emotion_items))
  ni <- nrow(obs_X)
  nrow_CovX <- ni * k
  all_person_X[i, 1:nrow_CovX] <- matrix(c(t(obs_X)), nrow = 1)
}


# Predict factor scores using OU factor model
factor_scores <- predict_eta(measurement_n, Lambda, Psi_list,
                             all_person_X, CovX_list)
factor_scores <- data.frame(factor_scores)
colnames(factor_scores) <- c('eta1', 'eta2') 
factor_scores <- factor_scores %>%
  mutate(user.id = rep(1:length(measurement_n), measurement_n))
dat_long <- dat_long %>%
  mutate(eta1 = factor_scores$eta1, eta2 = factor_scores$eta2)
# these factor scores, called eta1 and eta2 for the two latent factors,
#    could represent affective states of positive affect and negative affect

# Now that we've generated factor scores, we'll drop the observed long. outcomes
dat_long <- dat_long %>%
  dplyr::select(c(user.id, time, eta1, eta2))

# Combine the factor scores and the event information
dat_long <- dat_long %>% # NEW
  group_by(user.id) %>%
  mutate(start_time = time, stop_time = lead(time),
         obs_eta1_start = eta1, obs_eta2_start = eta2,
         obs_eta1_stop = lead(eta1), obs_eta2_stop = lead(eta2)) %>%
  ungroup()

dat_event <- read.csv(paste0(wd, 'data/sim_dat_event_setting', setting, '.csv'))
dat_cumulative <- left_join(dat_event, dat_long,
                            by = c('user.id' = 'user.id',
                                   'start_time' = 'start_time',
                                   'stop_time' = 'stop_time'))

# These new variables are defined as
# - obs_eta1_start = value of latent factor 1 at left endpoint of event interval
# - obs_eta1_stop = value of latent factor 1 at right endpoint of event interval
# - obs_eta2_start = value of latent factor 2 at left endpoint of event interval
# - obs_eta2_stop = value of latent factor 2 at right endpoint of event interval

