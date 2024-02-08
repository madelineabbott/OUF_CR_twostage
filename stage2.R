################################################################################
##                Fit cumulative risk model to simulated data                 ##
################################################################################

# Note: need to run stage 1 first

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
setting <- 1 # OU process has higher correlation
# setting <- 2 # OU process has lower correlation

# Set seed
g <- 1
set.seed(412+g)

# Set number of bootstrapped samples
B <- 200
# Set number of synthetic average values generated per event interval
R <- 50

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

# In stage 1, we fit the OUF model and predicted the factor scores
#. at measurement occasions

################################################################################
## Generate synthetic averages for latent factors and fit cumulative risk
## model to each bootstrapped dataset (ok if not interested in standard errors)
################################################################################

start_stage2_est <- Sys.time()

cat('-- Generating synthetic avg. values for entire dataset --\n')
b <- 0 # bootstrapped ID = 0 (not actually bootstrapped)
M <- 10 # number of multiple imputations
source(paste0(wd, 'aug_data.R'))

# Fit the weighted Poisson regression model
cat('-- Fitting weighted Poisson regression model for entire dataset --\n')
source(paste0(wd, 'fit_cr.R'))

end_stage2_est <- Sys.time()
results_rr$comp_sec <- as.numeric(end_stage2_est) - as.numeric(start_stage2_est)

# Save point estimates for this dataset
write.csv(x = results_rr,
          file = paste0(wd, 'stage2_results/nonbootstrapped_ests/results_setting',
                        setting, '_R', R, '.csv'),
          row.names = F)

################################################################################
## Bootstrap data and generate synthetic averages for latent factors
## And fit cumulative risk model to each bootstrapped dataset
################################################################################

M <- 2 # number of multiple imputations

start_stage2_inf <- Sys.time()

# Bootstrap original data and augment by generating R synthetic values for the
#  average of the latent factors across each event interval for all individuals
for (b in 1:B){
  cat('-- Generating synthetic avg. values for subsampled dataset #', b, '--\n')
  source(paste0(wd, 'aug_data_bs.R'))

  # Fit the weighted Poisson regression model
  cat('-- Fitting weighted Poisson regression model for #', b, '--\n')
  source(paste0(wd, 'fit_cr_bs.R'))
  # Save point estimates for this dataset
  write.csv(x = results_b2,
            file = paste0(wd, 'stage2_results/bootstrapped_ests/boot_results_setting',
                          setting, '_R', R, '_subsample', b, '.csv'),
            row.names = F)
  
}

end_stage2_inf <- Sys.time()

################################################################################
## Pool across bootstrapped results
################################################################################

# Calculate pooled point estimates and standard errors across the bootstraps
#  using the vonHippel approach
source(paste0(wd, 'combine_ests_cr_bs.R'))

final_results$comp_sec_rr <- as.numeric(end_stage2_est) - as.numeric(start_stage2_est)
final_results$comp_sec_vh <- as.numeric(end_stage2_inf) - as.numeric(start_stage2_inf)

# Save point estimates and standard errors
# est_vh and se_vh are results from bootstrap-based von Hippel approach 
# est_rr and se_rr are results from non-bootstrap Rubin's rule approach
write.csv(x = final_results,
          file = paste0(wd, 'stage2_results/estimates_setting', setting,
                        '_R', R, '.csv'),
          row.names = F)

# Save variance/covariance estimates (von Hippel only)
write.csv(x = VAR$var,
          file = paste0(wd, 'stage2_results/covariance_setting', setting,
                        '_R', R, '.csv'),
          row.names = F)


# ################################################################################
# ## Forest plot of point estimates and 95% confidence intervals
# ################################################################################
# 
# # Set values of true betas (cumulative risk model coefficients) used to simulate
# #  the data
# if (setting == 1){
#   final_results$true_values <- c(-2.4, -0.9,  1)
# } else if (setting == 2){
#   final_results$true_values <- c(-2.2, -0.6,  0.8)
# }
# 
# # True values are indicated with black circles in the plot while colored dots
# #  and error bars show point estimates and 95% confidence intervals
# ggplot(data = final_results) +
#   geom_errorbar(aes(ymin = est_vh - 1.96*se_vh,
#                     ymax = est_vh + 1.96*se_vh, x = variable,
#                     color = variable), width = 0.5) +
#   geom_point(aes(y = est_vh, x = variable, color = variable), size = 2) +
#   geom_point(aes(y = true_values, x = variable), size = 2, shape = 1) +
#   labs(x = 'Parameter', y = 'Estimate (95% CI)') +
#   scale_color_manual(values = pal) +
#   theme_bw() + theme(legend.position="none")


