################################################################################
# BOOTSTRAP DATA AND SAMPLE SYNTHETIC AVERAGE VALUES OF LATENT FACTORS FOR 
# FITTING THE CUMULATIVE RISK MODEL AND ESTIMATING STANDARD ERRORS
################################################################################

# This file bootstraps the simulated data and then augments each bootstrapped
#  subset with new synthetic values of the average latent process over each
#  event interval

################################################################################
# Step 1: Generate R synthetic values of the average of the latent process over
#  each event interval for this bootstrapped dataset
################################################################################

cat('   Step 1: generating sythetic values for avg. of latent process over event intervals \n')
#  Interval is (t, t + \Delta).  Ignore any information from Y(t) (event outcome)
all_etas <- matrix(nrow = 0, ncol = 8,
                   dimnames = list(NULL,
                                   c('i', 'ibs', 'm', 'r',
                                     'start_time', 'stop_time',
                                     'avg_eta1', 'avg_eta2')))



set.seed(214+b+g) # set new seed for boostrapped sample
# Bootstrap individual IDs
subsample_id <- sample(x = 1:N, size = N, replace =  TRUE)
# note that we don't use event/smoking info when generating sythetic average values

for (ibs in 1:length(subsample_id)){ # ibs = unique ID for bootstrapping
  i <- subsample_id[ibs]
  avg_etas_i <- data.frame(matrix(NA, nrow = 0, ncol = 6,
                                  dimnames = list(NULL, c('m', 'r',
                                                          'start_time', 'stop_time',
                                                          'avg_eta1', 'avg_eta2'))))
  
  cur_dat <- dat_cumulative %>%
      filter(user.id == i)
  ni <- nrow(cur_dat)
  eta1_anchors <- c(cur_dat$obs_eta1_start, cur_dat$obs_eta1_stop[ni])
  eta2_anchors <- c(cur_dat$obs_eta2_start, cur_dat$obs_eta2_stop[ni])
  
  # measurement times
  meas_times_i <- c(cur_dat$start_time, cur_dat$stop_time[ni])
  
  # Augment data with synthetic average values of latent factors, or etas
  # Von Hippel approach to standard error estimation uses multiple imputation
  #  within each bootstrapped sample
  # m indexes the imputation
  for (m in 1:M){
    #cat('  m =', m, '\n')
    # Assume etas at endpoints of intervals are known but avg eta is unknown
    # We directly sample the mean of the factor, given known values at endpoints
    cur_etas <- data.frame(time = meas_times_i, # times defining event intervals
                           eta1 = eta1_anchors, # factor 1 at endpoints
                           eta2 = eta2_anchors) # factor 2 at endpoints
    
    # Sample average eta R times
    aug_avg_etas <- sample_avg_etas_for_sims(obs_dat = cur_etas,
                                             theta = theta_boup,
                                             sigma = sigma_boup,
                                             R = R)
    aug_avg_etas <- aug_avg_etas %>%
      mutate(m = m) %>%
      dplyr::select(c('m', 'r',
                      'start_time', 'stop_time',
                      'avg_eta1', 'avg_eta2'))
    avg_etas_i <- rbind(avg_etas_i, aug_avg_etas)
    avg_etas_i <- avg_etas_i %>%
      arrange(r, start_time)
    
  } # end M loop
  
  # combine etas across all subjects
  all_etas <- rbind(all_etas, cbind(i, ibs, avg_etas_i))
}

cat('       Done with eta augmentation \n')

# restrict later analysis to all data after the start of the study ("post-quit")
all_etas <- all_etas %>%
  filter(start_time >= 0)

################################################################################
# Step 2: Incorporate cumulative cig counts
################################################################################

cat('Step 2: incorporating cumulative event counts \n')

all_etas <- all_etas %>%
  group_by(ibs, m, r) %>%
  mutate(smoke_int_id = 1:n(),
         avg_time = (start_time + stop_time)/2,
         smoking_interval_width = stop_time - start_time) %>%
  ungroup() %>%
  dplyr::select(c(ibs, i, m, r, smoke_int_id, avg_time,
                  smoking_interval_width, avg_eta1, avg_eta2))

# add event (cig info) to augmented etas for this boostrapped subset
true_cigs <- dat_cumulative %>% 
  dplyr::select(c(user.id, smoke_int_id, Y))
dat_cumulative_bs <- left_join(all_etas, true_cigs,
                               by = c('i' = 'user.id', 'smoke_int_id'))




 



