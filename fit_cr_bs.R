################################################################################
# FIT CUMULATIVE RISK MODEL
################################################################################

# PURPOSE: this file fits the cumulative risk model to the bootstrapped and
#  augmented data

# Set initialization approach: 'empirical' or 'uniform'
init_approach = 'uniform'

wt_list <- list()
coef_res <- matrix(nrow = 0, ncol = 4,
                   dimnames = list(NULL, c('iter', 'beta0', 'beta1', 'beta2')))

cat('   Starting estimation \n')

convg_crit <- FALSE

max_iter <- 150
for (iter in 1:max_iter){
  if (iter == 1){
    if (init_approach == 'empirical'){
      cat('   Initialization: intercept-only at empirical smoking rate \n')
      beta0 <- log(mean(dat_cumulative$Y/dat_cumulative$smoking_interval_width)) #true_betas[1]
      beta1 <- 0 
      beta2 <- 0 
      
      # step a: update weights
      dat_cumulative <- dat_cumulative %>%
        mutate(lambda_j = smoking_interval_width *
                 exp(beta0 + beta1*avg_eta1 + beta2*avg_eta2)) %>%
        mutate(PrY = exp(-lambda_j) * lambda_j ^ Y) %>%
        group_by(ibs, smoke_int_id) %>%
        mutate(total_PrY = sum(PrY)) %>%
        ungroup() %>%
        mutate(w = PrY / total_PrY) %>%
        mutate(w = w / sum(w) * N) # weight for each interval!
    }else if (init_approach == 'uniform'){
      cat('   Initialization: fit model with uniform weights \n')
      wtd_dat <- svydesign(ids = ~ ibs, weights = ~ 1,
                           data = dat_cumulative)
      mod0 <- svyglm(formula = Y ~ avg_eta1 + avg_eta2, family = quasipoisson,
                     offset = log(smoking_interval_width),
                     design = wtd_dat)

      beta0 <- coef(mod0)[which(names(coef(mod0)) == '(Intercept)')]
      beta1 <- coef(mod0)[which(names(coef(mod0)) == 'avg_eta1')]
      beta2 <- coef(mod0)[which(names(coef(mod0)) == 'avg_eta2')]

      # step a: update weights
      dat_cumulative <- dat_cumulative %>%
        mutate(lambda_j = smoking_interval_width *
                 exp(beta0 + beta1*avg_eta1 + beta2*avg_eta2)) %>%
        mutate(PrY = exp(-lambda_j) * lambda_j ^ Y) %>%
        group_by(ibs, smoke_int_id) %>%
        mutate(total_PrY = sum(PrY)) %>%
        ungroup() %>%
        mutate(w = PrY / total_PrY) %>%
        mutate(w = w / sum(w) * N) # weight for each interval!
    }

  }else{
    # fit weighted poisson model
    # for each m = 1, ..., M set of R imputed trajectories, fit model
    beta0_m <- rep(NA, M); beta1_m <- rep(NA, M); beta2_m <- rep(NA, M)
    beta0_se_m <- rep(NA, M); beta1_se_m <- rep(NA, M); beta2_se_m <- rep(NA, M)
    for (cur_m in 1:M){
      dat_cumulative_m <- dat_cumulative %>% filter(m == cur_m)
      wtd_dat <- svydesign(ids = ~ ibs, weights = ~ w,
                           data = dat_cumulative_m)
      mod0 <- svyglm(formula = Y ~ avg_eta1 + avg_eta2, family = quasipoisson,
                     offset = log(smoking_interval_width),
                     design = wtd_dat)
      beta0_m[cur_m] <- coef(mod0)[which(names(coef(mod0)) == '(Intercept)')]
      beta1_m[cur_m]  <- coef(mod0)[which(names(coef(mod0)) == 'avg_eta1')]
      beta2_m[cur_m]  <- coef(mod0)[which(names(coef(mod0)) == 'avg_eta2')]

      beta0_se_m[cur_m] <- sqrt(diag(vcov(mod0)))[1]
      beta1_se_m[cur_m] <- sqrt(diag(vcov(mod0)))[2]
      beta2_se_m[cur_m] <- sqrt(diag(vcov(mod0)))[3]

      # update weights
      dat_cumulative_m <- dat_cumulative_m %>%
        mutate(lambda_j = smoking_interval_width *
                 exp(beta0_m[cur_m] + beta1_m[cur_m]*avg_eta1 +
                       beta2_m[cur_m]*avg_eta2)) %>%
        mutate(PrY = exp(-lambda_j) * lambda_j ^ Y) %>%
        group_by(ibs, smoke_int_id) %>%
        mutate(total_PrY = sum(PrY)) %>%
        ungroup() %>%
        mutate(w = PrY / total_PrY) %>%
        mutate(w = w / sum(w) * N) # weight for each interval!

      # update full dataset with new weights
      dat_cumulative <- dat_cumulative %>% filter(m != cur_m)
      dat_cumulative <- rbind(dat_cumulative, dat_cumulative_m)
    }

    # update pooled coefficient estimates
    beta0 <- mean(beta0_m)
    beta1 <- mean(beta1_m)
    beta2 <- mean(beta2_m)
  }

  coef_res <- rbind(coef_res, c(iter, beta0, beta1, beta2))

  # check for convergence
  if (iter > 2){
    delta_ests <- max(abs(coef_res[iter, 2:4] - coef_res[iter-1, 2:4]))
    convg_crit <- ifelse(delta_ests < 1e-6, TRUE, FALSE)
    if(convg_crit == TRUE){
      cat('   Convergence in', iter, 'iterations: small change in estimates \n')
      break()
    }
  }

}
cat('   Done with point estimation \n')

# point estimates with SEs for this bootstrapped subsample of original data
results_b2 <- data.frame(B = b, M = 1:M,
                         beta0_m, beta1_m, beta2_m,
                        beta0_se_m, beta1_se_m, beta2_se_m)








