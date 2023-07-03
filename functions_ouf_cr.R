################################################################################
##  Functions used when fitting the dynamic Ornstein-Uhlenbeck factor model   ##
##              and the cumulative risk Poisson regression model              ##
################################################################################

library(Matrix)
library(mvtnorm)
library(reshape2)
library(dplyr)
library(ggplot2)
library(foreach)
library(cowplot)
library(msos)
library(lavaan)
library(semTools)
library(nlme)
library(parallel)

################################################################################
## FOR OU COVARIANCE MATRIX
################################################################################

# Calculate Kronecker sum
kronsum <- function(A, B){
  dimA <- nrow(A)
  dimB <- nrow(B)
  kronecker(A, diag(dimB)) + kronecker(diag(dimA), B)
}

# for marginal covariance (if eta(t_0) is unknown)
Gamma2 <- function(theta_ou, sigma_ou, t, s){ # t <= s
  if(t > s){
    return('error: t <= s is required')
  }
  theta_ks <- kronsum(theta_ou, theta_ou)
  theta_ks.inv <- solve(theta_ks)
  
  term1 <- theta_ks.inv
  term2 <- expm(t*theta_ks - kronsum(theta_ou*t, theta_ou*s))
  term3 <- matrix(sigma_ou %*% t(sigma_ou), ncol = 1)
  
  out <- matrix(term1 %*% term2 %*% term3,
                nrow = nrow(theta_ou), ncol = nrow(theta_ou))
  
  return(t(as.matrix(out,
                     row = nrow(theta_ou), ncol = nrow(theta_ou))))
}


# Calculate sigma using identifiability constraint that the OUP has a stationary
#  variance of 1
calc_sigma <- function(theta){
  p <- nrow(theta)
  theta_ks <- kronecker(theta, diag(p)) + kronecker(diag(p), theta)
  theta_ks.inv <- solve(theta_ks)
  sigma_vec <- c(diag(p))
  X <- theta_ks.inv[sigma_vec == 1, sigma_vec == 1]
  Y <- matrix(1, nrow = p, ncol = 1)
  sigma2_solutions <- c(solve(X, Y))
  return(sqrt(sigma2_solutions))
}


################################################################################
## FOR RE-SCALING OUP PARAMETERS TO SATISFY THE IDENTIFIABILITY CONSTRAINT
################################################################################

# Function for the covariance of the OU process (just the diagonal block)
construct_V <- function(theta, sigma) {
  kron_sum_theta = kronsum(theta, theta)  
  vec_sigmasq = as.vector(sigma %*% t(sigma))
  matrix(solve(kron_sum_theta, vec_sigmasq), nrow = nrow(theta))
}

# Project theta to identifiability space of theta*
update_theta <- function(theta, constants) {
  constants %*% theta %*% solve(constants)
}

# Project sigma to identifiability space of sigma*
update_sigma <- function(sigma, constants) {
  constants %*% sigma
}

# Calculate the constants that rescale the original OUP to one that satisfies
#   the identifiability constraint of constant stationary variance = 1
# Goal: want (c1, c2) s.t. cov2cor(Psi(theta, sigma)) = Psi(theta*, sigma*)
# Input: original OU parameters, guess of rescaling constants
# Output: rescaling constants (assuming > 0) and optimization convergence code
calc_constants <- function(guess_constants_vec, theta, sigma){
  log_constants_vec <- log(guess_constants_vec)
  obj_fun <- function(log_constants_vec, theta, sigma){
    constants <- diag(x = exp(log_constants_vec),
                      nrow = length(log_constants_vec))
    target_correlation <- cov2cor(construct_V(theta, sigma))
    theta_star <- update_theta(theta, constants)
    sigma_star <- update_sigma(sigma, constants)
    # cat('  ..... current sigma* is', diag(sigma_star), '\n')
    diag(sigma_star) <- pmax(diag(sigma_star), rep(0.001, nrow(sigma_star)))
    rescaled_covariance <- construct_V(theta_star, sigma_star)
    cur_diff <- sum(abs(target_correlation - rescaled_covariance))
    return(cur_diff)
  }
  if(nrow(theta)==1){
    opt <- optimize(interval = c(log(1e-4), log(99999)),
                    f = obj_fun, 
                    theta = theta, sigma = sigma)
    return(list(c_vec = exp(opt$minimum), code = 0))
  }else{
    opt <- optim(par = rep(0, length(log_constants_vec)), fn = obj_fun, theta = theta,
                 sigma = sigma)
    return(list(c_vec = exp(opt$par), code = opt$convergence))
  }
  
}
# example: calc_constants(c(1, 1), theta, sigma)


################################################################################
# ANALYTICAL GRADIENT (for measurement submodel block update)
################################################################################

# gradient w.r.t. Lambda
grad_lambda <- function(lambda_vec, Yi, Sigma_u, Theta, Psi){ #lambda_vec
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  lambda_vec[1] <- exp(lambda_vec[1])
  Lambda[nonzero==1] <- lambda_vec
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)
  j <- 1
  grad_vec <- rep(0, length(lambda_vec))
  
  for (i in 1:length(Lambda)){
    if (nonzero[i] == 1){
      J.k <- matrix(0, nrow = k, ncol = p)
      J.k[i] <- 1
      Ci <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(J.k))
      Ci <- Ci + t(Ci)
      inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
      term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
      term1 <- sum(diag(Sigma_star.inv %*% Ci))
      
      grad_vec[j] <- term1 + term2
      
      j <- j + 1
    }
  }
  ddlog <- rep(1, length(lambda_vec))
  ddlog[1] <- lambda_vec[1]
  return(-0.5*grad_vec*ddlog)
}

# gradient w.r.t. Sigma_u
grad_sigma_u <- function(R.sigma_u_vec, Yi, Lambda, Theta, Psi){
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp 
  
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)
  
  # grad_vec <- rep(0, length(R.sigma_u_vec)) # for non-diag Sigma_u
  grad_vec <- rep(0, k) # for diag Sigma_u
  
  # # for non-diag Sigma_u
  # j <- 1
  # for (i in 1:length(R.Sigma_u)){
  #   if (upper.tri(R.Sigma_u, diag = TRUE)[i] == 1){
  #     J.k <- matrix(0, nrow = k, ncol = k)
  #     J.k[i] <- 1
  #     ddR <- t(J.k) %*% R.Sigma_u + t(R.Sigma_u) %*% J.k
  #     Ci <- kronecker(ones_ni, ddR)
  # 
  #     term1 <- sum(diag(Sigma_star.inv %*% Ci))
  # 
  #     inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
  #     term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
  # 
  #     grad_vec[j] <- term1 + term2
  #     j <- j + 1
  #   }
  # }
  
  # for diagonal Sigma_u
  for (i in 1:k){
    J.k <- matrix(0, nrow = k, ncol = k)
    J.k[i,i] <- 1
    ddR <- t(J.k) %*% R.Sigma_u + t(R.Sigma_u) %*% J.k
    Ci <- kronecker(ones_ni, ddR)
    
    term1 <- sum(diag(Sigma_star.inv %*% Ci))
    
    inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
    term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
    
    grad_vec[i] <- term1 + term2
  }
  
  D <- matrix(0, nrow = k, ncol = k) # Note: switch 0 to 1 to estimate correlations
  diag(D) <- diag(R.Sigma_u)
  # D.vec <- D[upper.tri(D, diag = TRUE)] # for non-diag Sigma_u
  D.vec <- diag(D) # for diag Sigma_u
  return(-0.5 * grad_vec*D.vec)
}

# for Sigma_epsilon (called theta here)
grad_theta <- function(log.theta_vec, Yi, Lambda, Sigma_u, Psi){
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  
  Theta <- diag(exp(log.theta_vec)^2) # parameterized as log(sigma_e)
  
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)
  
  grad_vec <- rep(0, length(log.theta_vec))
  
  for (i in 1:length(log.theta_vec)){
    J.k <- matrix(0, nrow = k, ncol = k)
    J.k[i,i] <- 2*exp(log.theta_vec[i])
    # J.k[i,i] <- 2*(log.theta_vec[i]) 
    
    Ci <- kronecker(diag(ni), J.k)
    term1 <- sum(diag(Sigma_star.inv %*% Ci))
    inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
    term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
    grad_vec[i] <- term1 + term2
  }
  
  return(-0.5*grad_vec*exp(log.theta_vec))
  # return(-0.5*grad_vec) 
  
}

# gradient of the log-likelihood for i w.r.t. all meas submod parameters
lli_grad <- function(m_params, Yi, Psi_i){
  lambda_vec <- m_params[1:k]
  R.sigma_u_vec <- rep(0, choose(k, 2)+k)
  
  m = diag(k)
  upperm = upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  
  R.sigma_u_vec[diag(m)] <- m_params[(k+1):(2*k)]
  log.theta_vec <- m_params[(2*k+1):(3*k)]
  
  # set up Lambda
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero==1] <- lambda_vec
  cur_Lambda[1,1] <- exp(cur_Lambda[1,1])
  # set up Sigma_u
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  cur_Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp
  # set up Sigma_e (called Theta here)
  cur_Theta <- diag((exp(log.theta_vec))^2) # parameterized as log(sigma_e)
  # cur_Theta <- diag(((log.theta_vec))^2) 
  
  ni = length(Yi)/k
  Sigma_star_inv <- calc_Sigma_star_inv(ni, cur_Lambda, cur_Sigma_u,
                                        cur_Theta, Psi_i)
  
  c(c(grad_lambda_cpp_slow(k, ni, Yi, cur_Lambda, Psi_i, Sigma_star_inv,
                           nonzero)),
    c(grad_sigma_u_cpp_slow(k, ni, Yi, cur_Sigma_u, Sigma_star_inv)),
    c(grad_sigma_e_cpp_slow(k, ni, Yi, cur_Theta, Sigma_star_inv)))
}

# gradient for negative log likelihood w.r.t FA parameters
ll_grad <- function(m_params, all_person_Y_centered, Psi_list){
  # define non-zero indices of covariance matrices (for Sigma_u and Sigma_e)
  m <- diag(k)
  upperm <- upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  sigma_u_inds <- diag(m)+k
  sigma_e_inds <- (max(sigma_u_inds)+1):(max(sigma_u_inds)+k)
  new_grad_inds <- c(1:k, sigma_u_inds, sigma_e_inds)
  
  ########
  lambda_vec <- m_params[1:k]
  R.sigma_u_vec <- rep(0, choose(k, 2)+k)
  m = diag(k)
  upperm = upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  R.sigma_u_vec[diag(m)] <- m_params[(k+1):(2*k)]
  log.theta_vec <- m_params[(2*k+1):(3*k)]
  
  # set up Lambda
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero==1] <- lambda_vec
  cur_Lambda[1,1] <- exp(cur_Lambda[1,1])
  # set up Sigma_u
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  cur_Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp
  # set up Sigma_e
  cur_Sigma_e <- diag((exp(log.theta_vec))^2) # parameterized as log(sigma_e)
  ########
  
  # cat('using', n_cores, 'cores \n')
  # cl <- parallel::makeCluster(n_cores)
  # doParallel::registerDoParallel(cl)
  clusterExport(cl, c("measurement_times", "k", "lli_grad", "p", "nonzero"))
  grad_list <- foreach(i = 1:n_subj, .packages = 'mrabbott') %dopar% {
    ni <- length(measurement_times[i,!is.na(measurement_times[i,])])
    Y_i <- all_person_Y_centered[i, 1:(k*ni)]
    Psi_i <- Psi_list[[i]]
    Sigma_star_inv <- calc_Sigma_star_inv(ni, cur_Lambda, cur_Sigma_u,
                                          cur_Sigma_e, Psi_i)
    I_SYYt = diag(k*ni) - Sigma_star_inv %*% Y_i %*% t(Y_i)
    Sigma_term <-  I_SYYt %*% Sigma_star_inv
    fa_grads(k, ni, Y_i, cur_Lambda, cur_Sigma_u, cur_Sigma_e, Psi_i,
             Sigma_star_inv, nonzero, I_SYYt, Sigma_term)
  }
  # stopImplicitCluster()
  # parallel::stopCluster(cl)
  
  grad_mat <- matrix(unlist(grad_list), ncol = length(new_grad_inds),
                     byrow = TRUE)
  out <- -colSums(grad_mat)
  return(out)
}

################################################################################
#### Define FABOUP likelihood w.r.t FA parameters
################################################################################

# for single subject
FA_negllk_i <- function(m_params_vec, i){
  ni <- length(measurement_times[i,!is.na(measurement_times[i,])])
  # update Lambda
  lambda_vec <- m_params_vec[1:k]
  lambda_vec[1] <- exp(m_params_vec[1])
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero == 1] <- lambda_vec
  cur_Lambda_mat <- kronecker(diag(ni), cur_Lambda)
  # update Sigma_u
  log_sigma_u_vec <- m_params_vec[(k+1):(2*k)]
  # add lower bound to prevent est from getting stuck at -Inf
  log_sigma_u_vec <- pmax(log_sigma_u_vec, log(1e-4))
  cur_Sigma_u <- diag((exp(log_sigma_u_vec))^2)
  cur_Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), cur_Sigma_u)
  # update Theta
  log_theta_vec <- m_params_vec[(2*k+1):(3*k)]
  # add lower bound to prevent est from getting stuck at -Inf
  log_theta_vec <- pmax(log_theta_vec, log(1e-4))
  cur_Theta <- diag((exp(log_theta_vec))^2)
  cur_Theta_mat <- kronecker(diag(ni), cur_Theta)
  # update covariance matrix for Yi
  cur_Sigma_star <- cur_Lambda_mat %*% Psi_list[[i]]
  cur_Sigma_star <- cur_Sigma_star %*% t(cur_Lambda_mat)
  cur_Sigma_star <- cur_Sigma_star + cur_Sigma_u_mat + cur_Theta_mat
  cur_Sigma_star <- (cur_Sigma_star + t(cur_Sigma_star))/2
  
  # return negative log likelihood
  return(0.5 * (logdet(cur_Sigma_star)
                + tr(solve(cur_Sigma_star, empirical_CovY_list[[i]]))))
}

# negative log likelihood function (w.r.t FA parameters) for all subjects
FA_negllk <- function(m_params_vec, Psi_list){
  # likelihood
  out <- (sum(unlist(lapply(1:n_subj,
                            FUN = function(i){ FA_negllk_i(m_params_vec, i) }))))
  # gradient
  attr(out, 'gradient') <- ll_grad(m_params_vec, all_person_Y_centered,
                                   Psi_list)
  
  return(out)
}

################################################################################
#### Define FABOUP likelihood w.r.t OUP parameters
################################################################################

### 1. Define likelihood for Y w.r.t. to BOUP parameters
# NOTE: sigma is estimated on log scale
BOUP_neg_loglik_i <- function(params, c_vec, i) {
  # cat('person #:', i, '\n')
  theta_param_length = length(as.vector(theta_ou))
  theta_current = matrix(params[1:theta_param_length], nrow = nrow(theta_ou))
  sigma_current = diag(x = exp(params[(theta_param_length+1):length(params)]),
                       nrow = length(exp(params[(theta_param_length+1):length(params)])))
  
  # Rescale BOUP parameters to satisfy the identifiability constraint
  theta_current <- update_theta(theta_current,
                                diag(x = c_vec, nrow = length(c_vec)))
  sigma_current <- update_sigma(sigma_current,
                                diag(x = c_vec, nrow = length(c_vec)))
  
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec, ni = ni,
                                    times = meas_times_i)
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2
  
  # update Cov(Y_i)
  Lambda_mat <- kronecker(diag(ni), Lambda) %*% cur_Psi_i %*% kronecker(diag(ni), t(Lambda))
  Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), Sigma_u)
  Sigma_e_mat <- kronecker(diag(ni), Sigma_e)
  CovY_i <- Lambda_mat + Sigma_u_mat + Sigma_e_mat # Cov(Y_i) = Sigma^*_i
  CovY_i <- (CovY_i + t(CovY_i))/2
  
  # -loglik(Y)
  0.5 * (logdet(CovY_i) + tr(solve(CovY_i, empirical_CovY_list[[i]])))
}

BOUP_negloglik_n <- function(params, c_vec) {
  theta_param_length <- length(as.vector(theta_ou))
  theta_current <- matrix(params[1:theta_param_length], nrow = nrow(theta_ou))
  
  # check that theta corresponds to mean-reverting process
  theta_ok <- as.numeric(all(Re(eigen(theta_current)$values) > 0))
  
  if (theta_ok == 0){
    cat('warning: eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    sigma_current = diag(x = exp(params[(theta_param_length+1):length(params)]),
                         nrow = length(exp(params[(theta_param_length+1):length(params)])))
    current_c_fit <- calc_constants(c_vec, theta_current, sigma_current)
    current_c_vec <- current_c_fit$c_vec
    
    # return negative log likelihood
    return(sum(unlist(lapply(1:n_subj, FUN = function(i){
      BOUP_neg_loglik_i(params, current_c_vec, i) }))))
  }
}


################################################################################
#### Define FABOUP likelihood w.r.t ALL parameters
################################################################################

# for single subject
FABOUP_negllk_i <- function(m_params_vec, i){
  
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  # update Lambda
  lambda_vec <- m_params_vec[1:k]
  lambda_vec[1] <- exp(m_params_vec[1])
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero == 1] <- lambda_vec
  cur_Lambda_mat <- kronecker(diag(ni), cur_Lambda)
  # update Sigma_u
  log_sigma_u_vec <- m_params_vec[(k+1):(2*k)]
  cur_Sigma_u <- diag((exp(log_sigma_u_vec))^2)
  cur_Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), cur_Sigma_u)
  # update Sigma_e (called Theta here)
  log_theta_vec <- m_params_vec[(2*k+1):(3*k)]
  cur_Theta <- diag((exp(log_theta_vec))^2)
  cur_Theta_mat <- kronecker(diag(ni), cur_Theta)
  # update theta_OU
  theta_param_length = length(as.vector(theta_ou))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)], nrow = nrow(theta_ou))
  # update sigma_OU
  sigma_current = diag(x = exp(m_params_vec[(3*k+theta_param_length+1):length(m_params_vec)]),
                       nrow = p)
  # update OUP covariance matrix
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  ni = length(c(meas_times_i))
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec,
                                    ni = ni, times = c(meas_times_i))
  
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2
  
  # update covariance matrix for Yi
  cur_Sigma_star <- cur_Lambda_mat %*% cur_Psi_i
  cur_Sigma_star <- cur_Sigma_star %*% t(cur_Lambda_mat)
  cur_Sigma_star <- cur_Sigma_star + cur_Sigma_u_mat + cur_Theta_mat
  cur_Sigma_star <- (cur_Sigma_star + t(cur_Sigma_star))/2
  
  return(0.5 * (logdet(cur_Sigma_star)
                + tr(solve(cur_Sigma_star, empirical_CovY_list[[i]]))))
  
}

# for all i = 1, ..., N subjects
FABOUP_negloglik_n <- function(m_params_vec, calc_sigma = TRUE) {
  theta_param_length = length(as.vector(theta_ou))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)],
                         nrow = nrow(theta_ou))
  
  if (calc_sigma == TRUE){
    sigma_current = calc_sigma(theta_current)
    m_params_long <- c(m_params_vec, log(sigma_current))
  }else{
    m_params_long <- m_params_vec
  }
  
  # check mean reverting OUP
  theta_ok <- as.numeric(all(Re(eigen(theta_current)$values) > 0))
  
  if (theta_ok == 0){ 
    cat('warning: eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    sum(unlist(lapply(1:n_subj, FUN = function(i){ FABOUP_negllk_i(m_params_long, i) })))
  }
}


################################################################################
# FUNCTIONS FOR AUGMENTING LATENT PROCESS w/ SYNTHETIC VALUES OF ETA (not avg)
################################################################################

# conditional mean of OU process, eta(t) | \eta(s), \eta(u)
cond_mean <- function(eta_start, eta_stop, Psi, dt){
  v_1 <- Psi[1 : nrow(theta_boup), 1 : nrow(theta_boup)]
  psi_21 <- Psi[(nrow(theta_boup)+1) : (nrow(theta_boup)*2),
                1 : nrow(theta_boup)]
  psi_23 <- Psi[(nrow(theta_boup)+1) : (nrow(theta_boup)*2),
                (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  psi_13 <- Psi[1 : nrow(theta_boup),
                (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  psi_31 <- t(psi_13)
  v_3 <- Psi[(2*nrow(theta_boup)+1) : (nrow(theta_boup)*3),
             (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  
  term1 <- cbind(psi_21, psi_23)
  term2 <- solve(rbind(cbind(v_1, psi_13), cbind(psi_31, v_3)))
  term3 <- rbind(eta_start, eta_stop)
  
  # conditional mean for eta(t) | \eta(s), \eta(u)
  term1 %*% term2 %*% term3 
}

# conditional covariance of OU process, eta(t) | \eta(s), \eta(u)
cond_var <- function(eta_start, eta_stop, Psi){
  v_1 <- Psi[1 : nrow(theta_boup), 1 : nrow(theta_boup)]
  psi_21 <- Psi[(nrow(theta_boup)+1) : (nrow(theta_boup)*2),
                1 : nrow(theta_boup)]
  psi_23 <- Psi[(nrow(theta_boup)+1) : (nrow(theta_boup)*2),
                (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  psi_13 <- Psi[1 : nrow(theta_boup),
                (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  psi_31 <- t(psi_13)
  v_3 <- Psi[(2*nrow(theta_boup)+1) : (nrow(theta_boup)*3),
             (2*nrow(theta_boup)+1) : (nrow(theta_boup)*3)]
  psi_12 <- t(psi_21)
  psi_32 <- Psi[(2*nrow(theta_boup)+1) : (nrow(theta_boup)*3),
                (nrow(theta_boup)+1) : (nrow(theta_boup)*2)]
  v2 <- Psi[(nrow(theta_boup)+1) : (nrow(theta_boup)*2),
            (nrow(theta_boup)+1) : (nrow(theta_boup)*2)]
  
  term1 <- cbind(psi_21, psi_23)
  term2 <- solve(rbind(cbind(v_1, psi_13), cbind(psi_31, v_3)))
  term3 <- rbind(psi_12, psi_32)
  
  # conditional covariance for eta(t) | \eta(s), \eta(u)
  v2 - term1 %*% term2 %*% term3
}

# draw a value of eta between eta_start and eta_stop
draw_eta2 <- function(eta_start, eta_stop, Psi, dt){
  cur_cond_mean <- cond_mean(eta_start, eta_stop, Psi, dt)
  cur_cond_var <- cond_var(eta_start, eta_stop, Psi)
  cur_cond_var <- (cur_cond_var + t(cur_cond_var))/2
  mvtnorm::rmvnorm(1, mean = cur_cond_mean, sigma = cur_cond_var)
}

# eta augmentation with sheet (correlated)
augment_etas <- function(dat, theta_boup, sigma_boup){
  # dat should have a column for time and a column for each eta
  cur_obs <- dat[,c('time', paste0('eta', 1:nrow(theta_boup)))]
  obs_times <- sort(cur_obs$time)
  obs_etas <- cur_obs[,c(paste0('eta', 1:nrow(theta_boup)))]
  cur_max_gap = max(lag(cur_obs$time), na.rm = TRUE)
  while(cur_max_gap > max_gap) {
    for (j in 1:(length(obs_times)-1)){
      start_time = obs_times[j]
      stop_time = obs_times[j+1]
      eta_start = matrix(unlist(obs_etas[j,]), ncol = 1)
      eta_stop = matrix(unlist(obs_etas[j+1,]), ncol = 1)
      
      dt = stop_time - start_time
      
      # draw eta in the middle
      if(dt > max_gap){
        mid_time = dt/2 + start_time
        
        theta_t = t(theta_boup)
        sigma2_vec = matrix(sigma_boup %*% t(sigma_boup), ncol = 1)
        kron_sum_theta = kronsum(theta_boup, theta_boup)
        times = c(start_time, mid_time, stop_time)
        ni = length(times)
        Omega <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_boup, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec,
                                    ni = ni, times = times)
        
        Psi <- solve(Omega)
        new_eta <- draw_eta2(eta_start, eta_stop, Psi, dt)
        cur_obs <- rbind(cur_obs, c(mid_time, new_eta))
      }
      
    }
    # update observed eta times
    cur_obs <- cur_obs %>%
      arrange(time)
    obs_times <- sort(cur_obs$time)
    obs_etas <- cur_obs[,2:(nrow(theta_boup)+1)]
    cur_max_gap = max(cur_obs$time-lag(cur_obs$time), na.rm = TRUE)
  }
  return(cur_obs)
}


################################################################################
# FUNCTIONS FOR AUGMENTING LATENT PROCESS w/ SYNTHETIC VALUES OF AVG ETA 
################################################################################

# These functions use the exact distribution of avg_eta(s,t) | eta(t_L), eta(t_R)

# Variance function for average latent factors, eta, over interval from 0 to dt
var_ou_mean <- function(theta, sigma, dt){
  # calculate variance of the average of etas across interval of length dt
  V <- construct_V(theta, sigma)
  
  term1 <- solve((t(theta)*dt) %*% (t(theta)*dt))
  term2 <- expm(-t(theta)*dt)
  term3 <- diag(nrow(theta)) - t(theta)*dt
  
  term1b <- solve((theta*dt) %*% (theta*dt))
  term2b <- expm(-(theta)*dt)
  term3b <- diag(nrow(theta)) - theta*dt
  
  A <- matrix(V %*% (term1 %*% (term2 - term3)), nrow(theta), nrow(theta))
  B <- matrix(term1b %*% (term2b - term3b) %*% V, nrow(theta), nrow(theta))
  return(A + B)
}


# Covariance function for average eta over interval from 0 to dt;
# that is, [avg_eta(s,t), eta(t_L), eta(t_R)]^T
cov_ou_mean <- function(theta, sigma, time_L, time_s, time_t, time_R){
  # theta, sigma are ou parameters
  # t_left, t_right make endpoint of intervals where eta(tL), eta(tR) are known
  # s, t are endpoints of interval for avg eta
  V <- construct_V(theta, sigma)
  
  # Var(avg_eta(s,t))
  avg_eta_interval_width <- time_t - time_s
  var_avg_eta <- var_ou_mean(theta, sigma, avg_eta_interval_width)
  
  # Cov(eta(t_L), avg_eta(s, t))
  term1_L <- -1 / (time_t - time_s) * solve(theta)
  term2_L <- expm(-theta * (time_t - time_L)) - expm(-theta * (time_s - time_L))
  cov_etaL_avgeta <- term1_L %*% term2_L %*% V
  
  # Cov(avg_eta(s, t), eta(t_R))
  term1_R <- -1 / (time_t - time_s) * solve(t(theta))
  term2_R <- expm(-t(theta) * (time_R - time_s)) - expm(-t(theta) * (time_R - time_t))
  cov_avgeta_etaR <- V %*% term1_R %*% term2_R
  
  # Cov(eta(t_L), eta(t_R))
  cov_etaL_etaR <- expm(-theta * (time_R - time_L)) %*% V
  
  # Combine into single matrix: Cov(avgeta, eta_L, eta_R)
  row1 <- cbind(var_avg_eta, t(cov_etaL_avgeta), t(cov_avgeta_etaR))
  row2 <- cbind(cov_etaL_avgeta, V, cov_etaL_etaR)
  row3 <- cbind(cov_avgeta_etaR, t(cov_etaL_etaR), V)
  
  return(rbind(row1, row2, row3))
  
}

# Draw the average value of eta across the interval (s, t) given known values
#  of eta at the endpoints, time_L and time_R, where time_L <= s < t <= time_R
sample_avg_eta_j <- function(theta, sigma, time_L, time_s, time_t, time_R,
                             eta_L, eta_R, R){
  # etaL, etaR are both 2x1 matrices (or length-2 vecs) with current eta values
  # R is the number of draw wanted from avg_eta(s,t) | eta(t_L), eta(t_R)
  
  # covariance matrix of [avg_eta(s,t), eta(t_L), eta(t_R)]^T
  joint_covar <- cov_ou_mean(theta, sigma, time_L, time_s, time_t, time_R)
  # divide covariance matrix into blocks for calculating conditional covar
  block11 <- joint_covar[1:2, 1:2]
  block12 <- joint_covar[1:2, 3:6]
  block21 <- joint_covar[3:6, 1:2]
  block22 <- joint_covar[3:6, 3:6]
  
  # conditional mean of avg_eta | eta(t_L), eta(t_R)
  cond_mean <- block12 %*% solve(block22) %*% matrix(c(eta_L, eta_R), ncol = 1)
  cond_var <- as.matrix(block11 - block12 %*% solve(block22) %*% block21, 2, 2)
  # average to avoid small numerical errors if variance is tiny
  cond_var <- cond_var %*% t(cond_var) / 2
  
  # sample avg_eta(s,t) | eta(t_L), eta(t_R)
  avg_eta <- mvtnorm::rmvnorm(n = (R+1), mean = cond_mean, sigma = cond_var)
  res <- cbind(0:R, time_s, time_t, avg_eta)
  colnames(res) <- c('r', 'start_time', 'stop_time', 'avg_eta1', 'avg_eta2')
  return(res)
}


# This function is used for augmenting the data with synthetic avg eta values.
# Draw R+1 samples of the average value of eta acros the interval (s,t) given
#  known values of eta at the endpoints, time_L and time_R, where
#  time_L <= s < t <= time_R
# Use this function when analyzing the smoking cessation data
sample_avg_etas <- function(obs_dat, theta, sigma, R){
  # obs_dat is a dataframe with columns labeled 'time', 'eta1', 'eta2' that
  #  contains the information on the observed eta values
  # theta and sigma are OUP parameters
  # R is the number of sythetic trajectories (plus one with which is truth for simulations)
  avg_obs <- data.frame(matrix(NA, nrow = 0, ncol = 5,
                               dimnames = list(NULL,
                                               c('r', 'start_time', 'stop_time',
                                                 'avg_eta1', 'avg_eta2'))))
  # in avg_obs, start_time and stop_time correspond to time_s and time_t (ie the
  #  start and stop time for the augmented avg. eta interval)
  
  for (j in 1:nrow(obs_dat)){
    # define endpoints of interval of known etas
    time_L <- obs_dat$time_L[j]
    time_R <- obs_dat$time_R[j]
    time_s <- obs_dat$time_s[j]
    time_t <- obs_dat$time_t[j]
    eta_L <- c(unlist(obs_dat[j,c('eta1_start', 'eta2_start')]))
    eta_R <- c(unlist(obs_dat[j,c('eta1_stop', 'eta2_stop')]))
    # draw R values of the average of the augmented etas in this interval
    avg_etas2_direct <- sample_avg_eta_j(theta, sigma,
                                         time_L = time_L, time_s = time_s,
                                         time_t = time_t, time_R = time_R,
                                         eta_L = eta_L, eta_R = eta_R, R = R)
    # and update avg augmented eta dataset
    avg_obs <- rbind(avg_obs, avg_etas2_direct)
  }# end j loop
  return(avg_obs)
}

# This function does the same thing as the previous function but is just set
#  up slightly differently to handle simulated data.
# Use this function for augmenting the data with synthetic avg eta values.
sample_avg_etas_for_sims <- function(obs_dat, theta, sigma, R){
  # obs_dat is a dataframe with columns labeled 'time', 'eta1', 'eta2' that
  #  contains the information on the observed eta values
  # theta and sigma are OU process parameters
  # R is the number of sythetic trajectories (plus one with which is truth for simulations)
  avg_obs <- data.frame(matrix(NA, nrow = 0, ncol = 5,
                               dimnames = list(NULL,
                                               c('r', 'start_time', 'stop_time',
                                                 'avg_eta1', 'avg_eta2'))))
  # dat should have a column for time and a column for each eta
  obs_times <- obs_dat$time
  obs_etas <- obs_dat[,c('eta1', 'eta2')]
  
  for (j in 1:(length(obs_times)-1)){
    # define endpoints of interval of known etas
    time_L <- obs_times[j]
    time_R <- obs_times[j+1]
    eta_L <- c(unlist(obs_etas[j,]))
    eta_R <- c(unlist(obs_etas[j+1,]))
    # draw R values of the average of the augmented etas in this interval
    avg_etas2_direct <- sample_avg_eta_j(theta, sigma,
                                         time_L = time_L, time_s = time_L,
                                         time_t = time_R, time_R = time_R,
                                         eta_L = eta_L, eta_R = eta_R, R = R)
    # and update avg augmented eta dataset
    avg_obs <- rbind(avg_obs, avg_etas2_direct)
  }# end j loop
  return(avg_obs)
}

################################################################################
# FUNCTION FOR PLOTTING CROSS/AUTO-CORRELATION OF OU PROCESS OVER GAP TIMES
################################################################################

# Plot decay in correlation over time for OUP observed across different gap times
plot_OUP_cor <- function(theta_boup, sigma_boup, time_vec, plot_title){
  theta_boup_bonus <- theta_boup
  sigma_boup_bonus <- sigma_boup
  ## construct correlation matrix
  kron_sum_theta_bonus <- kronsum(theta_boup_bonus, theta_boup_bonus)
  sigma2_vec_bonus <- matrix(sigma_boup_bonus %*% t(sigma_boup_bonus), ncol = 1)
  ni <- length(time_vec)
  
  Omega_bonus <- calc_precision_cpp(kron_sum_theta = kron_sum_theta_bonus,
                                    theta = theta_boup_bonus, theta_t = t(theta_boup_bonus),
                                    sigma2_vec = sigma2_vec_bonus, ni = ni,
                                    times = time_vec)
  Psi_bonus <- cov2cor(solve(Omega_bonus))
  
  ## extract relevant blocks of correlation matrix
  # cor(eta_1(t), eta_1(s))
  cor_eta1t_eta1s_bonus <- data.frame(t_minus_s = time_vec,
                                      cor = Psi_bonus[1,seq(from = 1, to = ncol(Psi_bonus), by = 2)])
  
  # cor(eta_1(t), eta_2(s))
  cor_eta2t_eta2s_bonus <- data.frame(t_minus_s = time_vec,
                                      cor = Psi_bonus[2,seq(from = 2, to = ncol(Psi_bonus), by = 2)])
  
  # cor(eta_1(t), eta_2(s)), s <= t
  cor_eta1t_eta2s_bonus <- data.frame(t_minus_s = time_vec,
                                      cor = Psi_bonus[1,seq(from = 2, to = ncol(Psi_bonus), by = 2)])
  
  # cor(eta_1(a), eta_2(t)), s <= t
  cor_eta1s_eta2t_bonus <- data.frame(t_minus_s = time_vec,
                                      cor = Psi_bonus[2,seq(from = 1, to = ncol(Psi_bonus), by = 2)])
  
  ggplot() +
    geom_abline(slope = 0, intercept = 0, linetype = 'dotted') +
    geom_point(data = cor_eta1t_eta1s_bonus,
               aes(x = t_minus_s, y = cor, color = 'cor(eta1(t), eta1(s))'),
               size = 1, shape = 21) +
    geom_line(data = cor_eta1t_eta1s_bonus,
              aes(x = t_minus_s, y = cor, color = 'cor(eta1(t), eta1(s))')) +
    geom_point(data = cor_eta2t_eta2s_bonus,
               aes(x = t_minus_s, y = cor, color = 'cor(eta2(t), eta2(s))'),
               size = 1, shape = 22) +
    geom_line(data = cor_eta2t_eta2s_bonus,
              aes(x = t_minus_s, y = cor, color = 'cor(eta2(t), eta2(s))')) +
    geom_point(data = cor_eta1t_eta2s_bonus,
               aes(x = t_minus_s, y = cor, color = 'cor(eta1(t), eta2(s))'),
               size = 1, shape = 23) +
    geom_line(data = cor_eta1t_eta2s_bonus,
              aes(x = t_minus_s, y = cor, color = 'cor(eta1(t), eta2(s))')) +
    geom_point(data = cor_eta1s_eta2t_bonus,
               aes(x = t_minus_s, y = cor, color = 'cor(eta1(s), eta2(t))'),
               size = 1, shape = 24) +
    geom_line(data = cor_eta1s_eta2t_bonus,
              aes(x = t_minus_s, y = cor, color = 'cor(eta1(s), eta2(t))')) +
    theme_bw() +
    labs(x = '|t - s|, s <= t', y = 'correlation',
         subtitle = paste0(plot_title)) +
    scale_color_manual(values = c('purple', 'dodgerblue', 'orange', 'tomato'))
} 

################################################################################
# FUNCTION FOR PREDICTING FACTOR SCORES FROM LONGITUDINAL OU FACTOR MODEL
################################################################################

# Predict factor scores for 1 subject using dynamic (OU) factor model
predict_eta_i <- function(ni, Lambda, Psi_i, Y_i, Sigma_star_i){
  # ni = number of measurement occasions for subject i
  # Lambda = estimated loadings matrix
  # Psi_i = esitmated OU covariance matrix for subject i
  # Y_i = observed data for subject i (vector of length k x n_i)
  
  # and now predict factor scores
  Lambda_mat <- kronecker(diag(ni), Lambda)
  eta_hat <- Psi_i %*% t(Lambda_mat) %*% solve(Sigma_star_i) %*% Y_i
  
  # return factor scores
  return(matrix(eta_hat, ncol = p, byrow = T))
  
}

# Predict factor scores across all subjects
predict_eta <- function(measurement_n, Lambda, Psi_list, all_person_Y, CovY_list){
  factor_scores <- matrix(NA, nrow = 0, ncol = p)
  for (i in 1:length(measurement_n)){
    ni <- measurement_n[i]
    Psi_i <- Psi_list[[i]]
    Y_i <- all_person_Y[i, 1:(k*ni)]
    Sigma_star_i <- CovY_list[[i]]
    
    eta_i <- predict_eta_i(ni, Lambda, Psi_i, Y_i, Sigma_star_i)
    factor_scores <- rbind(factor_scores, eta_i)
  }
  
  return(factor_scores)
}
