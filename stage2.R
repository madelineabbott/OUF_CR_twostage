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
# Source Cpp code for calculating OU process precision matrix
Rcpp::sourceCpp(paste0(wd, "ou_precision_matrix.cpp"))

# Select set of true parameters under which data were simulated
setting <- 1 # OU process has high autocorrelation, low noise
# setting <- 2 # OU process has low autocorrelation, high noise

# Set number of bootstrapped subsamples
B <- 200
# Set number of multiple imputations
M <- 2
# Set number of synthetic average values generated per event interval
R <- 50

# Set seed
g <- 2
set.seed(412+g)
 
################################################################################
## Read in data
################################################################################

# dat_long contains info on the longitudinal outcomes
dat_long <- read.csv(paste0(wd, 'data/sim_dat_long_setting', setting, '.csv'))
# - user.id = ID for each individual
# - time = time of measurement occasion
# - active-scared = measured longitudinal outcomes

# dat_event contains info on the cumulative event outcome
dat_event <- read.csv(paste0(wd, 'data/sim_dat_event_setting', setting, '.csv'))
# - user.id = ID for each individual
# - start_time = time at left endpoint of event interval
# - avg_time = time at midpoint of event interval
# - stop_time = time at right endpoint of event interval
# - Y = cumulative number of events over this interval (e.g., # cig. smoked)
# - smoking_interval_width = width of event interval
# - smoke_int_id = within-individual ID for event intervals

# Set sample size (number of independent individuals)
N <- length(unique(dat_event$user.id))

################################################################################
## Define parameter values for OU factor model
################################################################################

# Define parameter values for the OU process--in the simulation studies, we use
#  the true OU process parameters, but these values could be replaced with
#  estimates (as in the data application)
if (setting == 1){
  # high autocorrelation, low noise
  theta_boup <- matrix(c(0.9, 0.5, 0.4, 1), 2, 2)
  sigma_boup <- diag(c(1.2, 1))
}else if (setting == 2){
  # low autocorrelation, high noise
  theta_boup <- matrix(c(1.8, 0.2, 0.4, 1.5), 2, 2)
  sigma_boup <- diag(c(4, 2))
}
p <- 2 # number of latent factors

# Define the other parameter matrices for the factor model
k <- 18 # number of longitudinal outcomes (# must match dat_long outcomes)
# Loadings matrix
Lambda <- matrix(c(0.821417185543206, 1.04458936956176,  0.876022383305528,
                   1.15900379444142,  0.810453790515167, 1.27250510014226, 
                   0.977819548205204, 1.14797926196996,  0.920847738140797,
                   rep(0, k),
                   1.18590850472731,  0.933236770307851, 1.24427973818947,
                   0.998974149217919, 1.16065757530552,  0.863549197572027,
                   1.01794611269844,  1.25405361079317,  1.05959968113476),
                 ncol = p)
# Covariance matrix for random intercepts
Sigma_u <- diag(rep(0, k))
# Covariance matrix for measurement error
Sigma_e <- diag(rep(0.1, k))

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
    dplyr::select(c('active', 'calm', 'determined', 'enthusiastic', 'grateful',
                    'happy', 'proud', 'joyful', 'attentive',
                    'angry', 'ashamed', 'disgusted', 'guilty', 'irritable',
                    'lonely', 'nervous', 'sad', 'scared'))
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

dat_cumulative <- left_join(dat_event, dat_long,
                            by = c('user.id' = 'user.id',
                                   'start_time' = 'start_time',
                                   'stop_time' = 'stop_time'))

# These new variables are defined as
# - obs_eta1_start = value of latent factor 1 at left endpoint of event interval
# - obs_eta1_stop = value of latent factor 1 at right endpoint of event interval
# - obs_eta2_start = value of latent factor 2 at left endpoint of event interval
# - obs_eta2_stop = value of latent factor 2 at right endpoint of event interval

################################################################################
## Bootstrap data and generate synthetic averages for latent factors
################################################################################

# Bootstrap original data and augment by generating R synthetic values for the
#  average of the latent factors across each event interval for all individuals
for (b in 1:B){
  cat('-- Generating synthetic avg. values for subsampled dataset #', b, '--\n')
  source(paste0(wd, 'aug_data_bs.R'))
  # save the boostrapped + augmented dataset
  write.csv(dat_cumulative_bs,
            file = paste0(wd, 'data/bootstrapped_data/sim_dat_setting', setting,
                          '_subsample', b, '.csv'))
}

################################################################################
## Fit cumulative risk model
################################################################################

# Fit weighted Poisson regression model to each bootstrapped dataset
for (b in 1:B){
  cat('-- Fit cumulative risk model to subsampled dataset #', b, '--\n')
  # Read in bootstrapped + augmented dataset
  dat_cumulative_bs <- read.csv(paste0(wd, 'data/bootstrapped_data/sim_dat_setting',
                                    setting, '_subsample', b, '.csv'))
  # Fit the weighted Poisson regression model
  source(paste0(wd, 'fit_cr_bs.R'))
  # Save point estimates for this dataset
  write.csv(x = results_b2,
            file = paste0(wd, 'data/bootstrapped_ests/boot_results_setting',
                          setting, '_subsample', b, '.csv'),
            row.names = F)
}

# Calculate pooled point estimates and standard errors across the bootstraps
#  using the vonHippel approach
source(paste0(wd, 'combine_ests_cr_bs.R'))

# Save point estimates and standard errors
write.csv(x = final_results,
          file = paste0(wd, 'data/estimates_setting', setting, '.csv'),
          row.names = F)

# Save variance/covariance estimates
write.csv(x = VAR$var,
          file = paste0(wd, 'data/covariance_setting', setting, '.csv'),
          row.names = F)


################################################################################
## Forest plot of point estimates and 95% confidence intervals
################################################################################

# Set values of true betas (cumulative risk model coefficients) used to simulate
#  the data
if (setting == 1){
  final_results$true_values <- c(-2.4, -0.9,  1)
} else if (setting == 2){
  final_results$true_values <- c(-2.2, -0.6,  0.8)
}

# True values are indicated with black circles in the plot while colored dots
#  and error bars show point estimates and 95% confidence intervals
ggplot(data = final_results) +
  geom_errorbar(aes(ymin = est - 1.96*se, ymax = est + 1.96*se, x = variable,
                    color = variable), width = 0.5) +
  geom_point(aes(y = est, x = variable, color = variable), size = 2) +
  geom_point(aes(y = true_values, x = variable), size = 2, shape = 1) +
  labs(x = 'Parameter', y = 'Estimate (95% CI)') +
  scale_color_manual(values = pal) +
  theme_bw() + theme(legend.position="none")


