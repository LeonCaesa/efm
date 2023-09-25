source('./R/utils.R')
source('./R/efm.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

set.seed(1)
# [fixed simulation configuration ]
adam_control <- adam.control(max_epoch = 5, batch_size = 128,
                             step_size = 1, rho = 0, abs_tol = 1e-6,
                             beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control <- sample.control(sample_size = 500, eval_size = 500)
fagqem_control <- list(maxit = 5, epsilon = 1e-6, trace = FALSE)
n_gq <- 15

q <- 3; dispersion_star <- 2
factor_weights <- 1
L_prior <- list(mean = rep(0, q),
                precision = rep(1, q))
V_prior <- list(mean = rep(0, q),
                sigma = rep(0.5, q) )


factor_family <- quasipoisson()
n <- 500; d_list <- seq(50, 100, by =10)
fagqem_control$maxit <- adam_control$max_epoch * as.integer(n / adam_control$batch_size)

for (d in d_list){
  center_star <- rep(0, d)
  truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                        weights = factor_weights, phi = dispersion_star)

  # [intialize script]
  init <- list(center = truth$center + rnorm(d, 0, 1),
               dispersion = rnorm(d, 0, truth$phi/3),
               Vt = svd(truth$X, nu = q, nv = q)$v)
  names_algo <- c('lapl', 'sml', 'ps')

  for (algo_ in names_algo){
    time_efm<- system.time(efm_result <- efm(truth$X/factor_weights, factor_family,
                                             rank = q, weights = factor_weights,
                                             start = init, algo = algo_,  adam_control = adam_control,
                                             sample_control = sample_control, eval_likeli = TRUE,
                                             lambda_prior = L_prior))
    efm_result$exe_time <- time_efm[3]
  } # end of algo loop

  time_fagqem <- system.time(
    efm_fagqem <- fa_gqem(truth$X, q, n_gq, family = factor_family,
                          start = init, control = fagqem_control, eval_likeli = TRUE,
                          eval_size = sample_control$eval_size)
  )
  efm_fagqem$exe_time <- time_fagqem[3]
  break

}


true_cov <- marg_var(truth$center, truth$V, family = factor_family, ngq = 15, dispersion = truth$phi)
efm_esti <- marg_var(efm_result$center, efm_identify(efm_result$V), family = factor_family, ngq = 15, dispersion = efm_result$dispersion)
fagqem_esti <- marg_var(efm_fagqem$center, efm_fagqem$V, family = factor_family, ngq = 15, dispersion = efm_fagqem$dispersion)
naive_esti <- cov(truth$X)

mse_naive = mean((naive_esti - true_cov)^2)
mse_efm = mean((efm_esti - true_cov)^2)
mse_fagqem = mean( (fagqem_esti - true_cov)^2 )
c(mse_naive, mse_efm, mse_fagqem)



