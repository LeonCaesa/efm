setwd('/projectnb/dmfgrp/efm/')
source('./R/utils.R')
source('./R/efm.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

set.seed(1)
# [fixed simulation configuration ]
save_dir = '/projectnb/dmfgrp/efm/CovResult/'

adam_control <- adam.control(max_epoch = 20, batch_size = 128,
                             step_size = 0.5, rho = 0, abs_tol = 1e-6,
                             beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control <- sample.control(sample_size = 500, eval_size = 500)
fagqem_control <- list(maxit = 5, epsilon = 1e-6, trace = FALSE)
n_gq <- 15

q <- 3; dispersion_star <- 2
factor_weights <- 1

# [2008 Jianqing.Fan Factor Model Params]
L_prior <- list(mean = c(0.023558, 0.012989, 0.020714),
                precision =  1/ c(1.2507, 0.31564, 0.19303))
V_prior <- list(mean = c(0.78282, 0.51803, 0.41003),
                sigma = matrix( c(c(0.029145, 0.023873, 0.010184),
                                c(0.023873, 0.053951, -0.006967),
                                c(0.010184, -0.006967, 0.086856)), nrow = 3))

phi_prior<-list(alpah = 4.0713, beta= 0.1623)




#factor_family <- quasipoisson()
factor_family <- negative.binomial(dispersion_star)
n <- 756; d_list <- seq(16, 1000, by = 100)
fagqem_control$maxit <- adam_control$max_epoch * as.integer(n / adam_control$batch_size)




for (d in d_list){
  center_star <- rep(0, d)
  dispersion_star <- rgamma(d, shape = phi_prior$alpah, scale =phi_prior$beta)^2


  truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                        weights = factor_weights, phi = dispersion_star)

  save_name <- paste(paste(save_dir, 'truth', d, n, sep = '_'), '.RData', sep = '')
  save(truth, file = save_name)


  # [intialize script]
  init <- NULL
  # init <- list(center = truth$center + rnorm(d, 0, 1),
  #              dispersion = rnorm(d, 0, truth$phi/3),
  #              Vt = svd(truth$X, nu = q, nv = q)$v)
  # save_name <- paste(paste(save_dir, 'init', d, n, sep = '_'), '.RData', sep = '')
  # save(init, file = save_name)

  time_fagqem <- system.time(
    efm_fagqem <- fa_gqem(truth$X, q, n_gq, family = factor_family, lambda_prior= L_prior,
                          start = init, control = fagqem_control, eval_likeli = TRUE,
                          eval_size = sample_control$eval_size)
  )
  efm_fagqem$exe_time <- time_fagqem[3]
  save_name <- paste(paste(save_dir, 'fageqm', d, n, sep = '_'), '.RData', sep = '')
  save(efm_fagqem, file = save_name)


  names_algo <- c('lapl', 'sml', 'ps')

  for (algo_ in names_algo){
    time_efm<- system.time(efm_result <- efm(truth$X/factor_weights, factor_family,
                                             rank = q, weights = factor_weights,
                                             start = init, algo = algo_,  adam_control = adam_control,
                                             sample_control = sample_control, eval_likeli = TRUE,
                                             lambda_prior = L_prior))
    efm_result$exe_time <- time_efm[3]

    save_name <- paste(paste(save_dir, 'efm', algo_, d, n, sep = '_'), '.RData', sep = '')
    save(efm_result, file = save_name)
  } # end of algo loop



}
#
#
# true_cov <- marg_var(truth$center, truth$V, family = factor_family, ngq = 15, dispersion = truth$phi)
# efm_esti <- marg_var(efm_result$center, efm_identify(efm_result$V), family = factor_family, ngq = 15, dispersion = efm_result$dispersion)
# fagqem_esti <- marg_var(efm_fagqem$center, efm_fagqem$V, family = factor_family, ngq = 15, dispersion = efm_fagqem$dispersion)
# naive_esti <- cov(truth$X)
#
# mse_naive = mean((naive_esti - true_cov)^2)
# mse_efm = mean((efm_esti - true_cov)^2)
# mse_fagqem = mean( (fagqem_esti - true_cov)^2 )
# c(mse_naive, mse_efm, mse_fagqem)
#
#
#
