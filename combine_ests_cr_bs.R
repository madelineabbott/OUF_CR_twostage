################################################################################
##         Pool results across the fitted cumulative risk models              ##
################################################################################

# This file combines estimates from the cumulative risk model fit to the
#  bootstrapped and augmented datasets following the von Hippel approach.

################################################################################
## Read in results from fitted models (bootstrapped)
################################################################################

results_allB <- matrix(nrow = 0, ncol = 5,
                       dimnames = list(NULL, c('B', 'M', 'beta0', 'beta1',
                                               'beta2')))

for (b in 1:B){
  results_b <- read.csv(paste0(wd, 'stage2_results/bootstrapped_ests/boot_results_setting',
                               setting, '_R', R, '_subsample', b, '.csv'))
  results_b <- results_b %>%
    dplyr::select(c('B', 'M', paste0('beta', 0:2, '_m')))
  results_allB <- rbind(results_allB, results_b)
}

colnames(results_allB) <- c('B', 'M', paste0('beta', 0:2))

################################################################################
## Point estimates
################################################################################

cat('   Calculating pooled point estimates \n')

beta0_BM <- mean(results_allB$beta0)
beta1_BM <- mean(results_allB$beta1)
beta2_BM <- mean(results_allB$beta2)

################################################################################
## Variance/covariance estimates (using von Hippel approach)
################################################################################

cat('   Calculating pooled variance estimates \n')

# for calculating mean sum of squares (w/in & b/n bootstraps)
calc_mss <- function(coef_names, results){
  ests <- results %>%
    dplyr::select(all_of(coef_names)) %>%
    as.matrix()
  
  preds <- matrix(NA, nrow = nrow(ests), ncol = ncol(ests),
                  dimnames = list(NULL, colnames(ests)))
  for (j in coef_names){
    # fit linear regression models
    fit_lm <- lm(results[,which(colnames(results) == j)] ~ factor(B),
                 data = results)
    preds[,which(colnames(preds) == j)] <- fitted(fit_lm)
  }
  
  # calculate residuals (between & within)
  E_matrix <- preds - matrix(rep(colMeans(ests), nrow(ests)),
                             nrow = nrow(ests), byrow = T)
  R_matrix <- preds - ests
  
  # mean sum of squares between bootstraps
  SSE <- t(E_matrix) %*% E_matrix
  MSB <- SSE / (n - b)
  
  # mean sum of squares within bootstraps
  SSR <- t(R_matrix) %*% R_matrix
  MSW <- SSR / (b - 1)
  
  # return MSWithin & MSBetween
  return(list(MSB = MSB,
              MSW = MSW))
}


# calculate variance using von Hippel approach
calc_var <- function(MSB, MSW){
  var_bs <- (B + 1) / (B*M) * MSB - MSW / M
  # cat('Is MSB - MSW < 0? ', sum(diag(MSB - MSW) < 0)>0, '\n')  
  
  df_num <- ((B+1)/(B*M) * MSB - MSW)^2
  df_den <- ((B+1)/(B*M))^2 * MSB^2 / (B-1) + MSW^2 / (B*M^2*(M-1))
  df <- matrix(df_num / df_den, nrow = nrow(df_den), ncol = ncol(df_den))
  
  # return variance estimate and Satterthwaite's df (these dfs are not used)
  list(var = var_bs, df = df)
}


# calculate covariance estimates
n <- nrow(results_allB)
b <- length(unique(results_allB$B)) # number of groups (bootstraps)

coef_names <- paste0('beta', 0:2)
MSS <- calc_mss(coef_names, results = results_allB)
VAR <- calc_var(MSB = MSS$MSB, MSW = MSS$MSW)

################################################################################
## Point and variance/covariance estimates (using Rubin's Rule)
################################################################################

# Read in results from fitted models
results_RR <- read.csv(paste0(wd, 'stage2_results/nonbootstrapped_ests/results_setting',
                              setting, '_R', R, '.csv'))


################################################################################
# Save results
################################################################################

final_results <- data.frame(variable = paste0('beta', 0:2),
                            est_vh = c(beta0_BM, beta1_BM, beta2_BM),
                            se_vh = sqrt(diag(VAR$var)), # von Hippel
                            est_rr = results_RR$est_rr,
                            se_rr = results_RR$se_rr) # Rubin's rule



