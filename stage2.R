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
# setting <- 1 # high autocorrelation, low noise
setting <- 2 # low autocorrelation, high noise

# Set number of bootstrapped subsamples
B <- 200
# Set number of multiple imputations
M <- 2
# Set number of synthetic average values generated per event interval
R <- 50

# Set seed
g <- 0

################################################################################
## Read in data
################################################################################
 
dat_cumulative <- read.csv(paste0(wd, 'data/sim_dat_setting', setting, '.csv'))

# variables in this dataset are defined as:
# - user.id = ID for each individual
# - start_time = left endpoint of event interval
# - avg_time = time at midpoint of event interval
# - stop_time = right endpoint of event interval
# - Y = cumulative # of events per interval (e.g., num. cigs. smoked)
# - smoking_interval_width = width of event interval
# - smoke_int_id = within-individual ID for event intervals
# - obs_eta1_start = value of latent factor 1 at left endpoint of event interval
# - obs_eta1_stop = value of latent factor 1 at right endpoint of event interval
# - obs_eta2_start = value of latent factor 2 at left endpoint of event interval
# - obs_eta2_stop = value of latent factor 2 at right endpoint of event interval

# Set sample size (number of individuals)
N <- length(unique(dat_cumulative$user.id))

################################################################################
## Bootstrap data and generate sythetic averages for latent factors
################################################################################

# Define parameter values for the OU process--in the simulation studies, we use
#  the true OU process parameters, but these values could be replaced with
#  estimates
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

# Bootstrap original data and augment by generating R sythetic values for the
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

# Fit weighted Poisson regression model to each boostrapped dataset
for (b in 1:B){
  cat('-- Fit cumulative risk model to subsampled dataset #', b, '--\n')
  # Read in boostrapped + augmented dataset
  dat_cumulative <- read.csv(paste0(wd, 'data/bootstrapped_data/sim_dat_setting',
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

# Set values of true betas (cumul. risk coeffs) to add to plot
if (setting == 1){
  final_results$true_values <- c(-2.4, -0.9,  1)
}else if (setting == 2){
  final_results$true_values <- c(-2.2, -0.6,  0.8)
  
}

# True values are indicated with black circles in the plot while colored dots and
#  error bars show estimates
ggplot(data = final_results) +
  geom_errorbar(aes(ymin = est - 1.96*se, ymax = est + 1.96*se, x = variable,
                    color = variable), width = 0.5) +
  geom_point(aes(y = est, x = variable, color = variable), size = 2) +
  geom_point(aes(y = true_values, x = variable), size = 2, shape = 1) +
  labs(x = 'Parameter', y = 'Estimate (95% CI)') +
  scale_color_manual(values = pal) +
  theme_bw() + theme(legend.position="none")


