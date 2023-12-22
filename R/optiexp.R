setwd('/projectnb/dmfgrp/efm/')
source('./R/utils.R')
source('./R/efm.R')
library(mvtnorm)
library(matrixStats)
library(MASS)



#Get aruments from the command line
argv <- commandArgs(TRUE)

# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
  family_idx <- as.numeric( argv[1] )
  sample_idx <- as.numeric( argv[2] )
  algo_idx <- as.numeric( argv[3] )
  d <- as.numeric( argv[4] )
  q <- as.numeric( argv[5] )
}
# family_idx = 3; sample_idx = 3; algo_idx = 1; d = 50; q = 2
# Print the values
paste("family index:", family_idx)
paste("sample_idx:", sample_idx)
paste("algo_idx:", algo_idx)
paste("d:", d)
paste("q:", q)


# [experiment specific]
sample_list = c(50, 300, 500)
algo_list = c('ps', 'sml', 'lapl', 'em')

# [family specifics]
dispersion_list = c(1, 1, 20)
family_list = c()
family_list[[1]] = poisson('log')
family_list[[2]] = Gamma('log')
family_list[[3]] =  negative.binomial(20)
name_list = c('poisson', 'Gamma', 'negbinom')
# [setup family]
factor_family = family_list[[family_idx]]
dispersion = dispersion_list[family_idx]

# [generate data]
set.seed(2)
n = 500
factor_weights = 1

L_prior <- list(mean = rep(0, q),
                precision =  rep(1,q))

V_prior <- list(mean = rep(0.5, q),
                sigma = rep(1, q))
center_star <- rep(0, d)
dispersion_star <- dispersion_list[family_idx]

truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                      weights = factor_weights, phi = dispersion_star)


# [optimization setup]
adam_control <- list(max_epoch = 25, batch_size = 128,
                     rho = 0, abs_tol =1e-6,
                     step_size = 0.5,
                     beta1 = 0.9, beta2 = 0.999,
                     epsilon = 1e-8)

sample_control <- list(sample_size = sample_list[sample_idx],
                       eval_size = 1500)

init <- list(center = truth$center + rnorm(d, 0, 1),
             dispersion = 1,
             Vt = svd(truth$X, nu = q, nv = q)$v)

save_name = paste( '/projectnb/dmfgrp/efm/OptiResult1221/',
                   paste( algo_list[algo_idx],
                     name_list[family_idx], paste('s', sample_list[sample_idx], sep =''),
                          paste('d', d, sep = ''),
                          paste('q', q, sep = ''),
                          paste('T', adam_control$max_epoch, sep= ''), sep = '_'),
                   '.RData', sep ='')


if (!file.exists(save_name)){
    set.seed(1)
    start = Sys.time()
    efm_result <- efm(x = truth$X/factor_weights, lambda_prior = L_prior,
                       factor_family = factor_family,
                       rank = q, weights = factor_weights,
                       algo = algo_list[algo_idx],
                       start = init,
                       sample_control = sample_control,
                       adam_control = adam_control,
                       eval_likeli = TRUE)

    end = Sys.time()
    efm_time = as.numeric(difftime(end, start, units = 's')) - efm_result$eval_time
    efm_result$efm_time = efm_time
    save(efm_result, file = save_name)
}




