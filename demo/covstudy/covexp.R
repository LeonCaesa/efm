library(efm)
#setwd('/projectnb/dmfgrp/efm/')
#source('../R/utils.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

set.seed(1)
# [2008 Jianqing.Fan Factor Model Params]
q <- 3 ;n <- 756
L_prior <- list(mean = c(0.023558, 0.012989, 0.020714),
                precision =  1/ c(1.2507, 0.31564, 0.19303))

V_prior <- list(mean = c(0.78282, 0.51803, 0.41003),
                sigma = matrix( c(c(0.029145, 0.023873, 0.010184),
                                  c(0.023873, 0.053951, -0.006967),
                                  c(0.010184, -0.006967, 0.086856)), nrow = 3))
phi_prior<-list(alpah = 4.0713, beta= 0.1623)

# [Optimization configuration]
n_gq <- 15
adam_control <- adam.control(max_epoch = 20, batch_size = 128,
                             step_size = 0.2, rho = 0, abs_tol = 1e-8,
                             beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8)
sample_control <- sample.control(sample_size = 500, eval_size = 500)
fagqem_control <- list(maxit = adam_control$max_epoch * as.integer(n / adam_control$batch_size),
                       epsilon = 1e-8, trace = FALSE)



# [Bash script configuration]
# Check if the command line is not empty and convert values to numerical values
argv <- commandArgs(TRUE)
if (length(argv) > 0){
  exp_idx <- argv[1]
  n_repeats <- as.numeric( argv[2] )
  d <- as.numeric( argv[3] )
}
# exp_idx = 2; n_repeats = 1; d = 16
print(c(exp_idx, n_repeats, d))
set.seed(n_repeats)


switch(exp_idx,
       "1" = {
         dispersion_star <- 20; factor_weights <- 1
         factor_family <- negative.binomial(dispersion_star)},
       "2" = {
         dispersion_star <- NULL; factor_weights <- 1
         factor_family <- quasipoisson()},
       "3" = {
         dispersion_star <- 1; factor_weights <- rpois(n*d, 20) + 1; dim(factor_weights) = c(n,d)
         factor_family <- binomial()},
       "4" = {
         dispersion_star <- 1; factor_weights <- 1
         factor_family <- poisson()},
       stop(paste0("No handler for ", exp_id))
)


#save_dir = paste('/projectnb/dmfgrp/efm/CovResult1209/', factor_family$family, '/', sep = '')
#dir.create(file.path(save_dir), showWarnings = FALSE)



# [main script below]
center_star <- rep(0, d)
if (check_DispersionUpdate(factor_family)){
  dispersion_star <- rgamma(d, shape = phi_prior$alpah, scale =phi_prior$beta)^2 + 1}

truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                      weights = factor_weights, phi = dispersion_star)

save_name <- paste(paste(save_dir, 'truth', d, n_repeats, sep = '_'), '.RData', sep = '')
save(truth, file = save_name)


# [initialize script]
init <- NULL
# init <- list(center = truth$center + rnorm(d, 0, 1),
#              dispersion = rnorm(d, 0, truth$phi/3),
#              Vt = svd(truth$X, nu = q, nv = q)$v)
# save_name <- paste(paste(save_dir, 'init', d, n, sep = '_'), '.RData', sep = '')
# save(init, file = save_name)

time_fagqem <- system.time(
  efm_fagqem <- fa_gqem(X = truth$X/factor_weights, q, n_gq, weights = factor_weights,
                        family = factor_family, lambda_prior= L_prior,
                        start = init, control = fagqem_control, eval_likeli = TRUE,
                        eval_size = sample_control$eval_size)
)
efm_fagqem$exe_time <- time_fagqem[3]
save_name <- paste(paste(save_dir, 'fageqm', d, n_repeats, sep = '_'), '.RData', sep = '')
save(efm_fagqem, file = save_name)


# names_algo <- c('lapl', 'sml', 'ps')
#
# for (algo_ in names_algo){
#   print(c(d, algo_))
#   time_efm<- system.time(efm_result <- efm(truth$X/factor_weights, factor_family,
#                                            rank = q, weights = factor_weights,
#                                            start = init, algo = algo_,  adam_control = adam_control,
#                                            sample_control = sample_control, eval_likeli = TRUE,
#                                            lambda_prior = L_prior))
#
#   efm_result$exe_time <- time_efm[3]
#
#   save_name <- paste(paste(save_dir, 'efm', algo_, d, n_repeats, sep = '_'), '.RData', sep = '')
#   save(efm_result, file = save_name)
# } # end of algo loop
