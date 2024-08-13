library(efm)
#setwd('/projectnb/dmfgrp/efm/')
#source('../R/utils.R')
library(mvtnorm)
library(matrixStats)
library(MASS)



#Get aruments from the command line
argv <- commandArgs(TRUE)

# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
  family_idx <- as.numeric( argv[1] )
  #sample_idx <- as.numeric( argv[2] )
  algo_idx <- as.numeric( argv[2] )
  d <- as.numeric( argv[3] )
  q <- as.numeric( argv[4] )
}
# family_idx = 1; algo_idx = 4; d = 512; q = 2
# Print the values
paste("family index:", family_idx)
#paste("sample_idx:", sample_idx)
paste("algo_idx:", algo_idx)
paste("d:", d)
paste("q:", q)


# [experiment specific]
sample_list = c(50, 300, 500)
algo_list = c('ps', 'sml', 'lapl', 'em')
name_list = c('poisson', 'binomial', 'negbinom')
dispersion_list = c(1, 1, 20)

# [family specifics]
family_list = c()
family_list[[1]] = poisson('log')
family_list[[2]] = binomial()
family_list[[3]] =  negative.binomial(20)
# [setup family]
factor_family = family_list[[family_idx]]
dispersion = dispersion_list[family_idx]



# [generate data]
set.seed(2)
n = 500
if (family_idx == 2){
  factor_weights <- rpois(n*d, 20) + 1
  dim(factor_weights) = c(n,d)
}else{
  factor_weights <- 1
}
# factor_weights <- 1


L_prior <- list(mean = rep(0, q),
                precision =  rep(1,q))

V_prior <- list(mean = rep(0, q),
                sigma = rep(0.3, q))
                #sigma = rep(1, q))
center_star <- rep(2, d)
dispersion_star <- dispersion_list[family_idx]

truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                      weights = factor_weights, phi = dispersion_star)
plot(density(c(truth$X)))

# [optimization setup]
adam_control <- list(max_epoch = 25, batch_size = 128,
                     rho = 0, abs_tol =1e-6,
                     step_size = 0.03,
                     beta1 = 0.9, beta2 = 0.999,
                     epsilon = 1e-8)



em_control <- list(maxit = 25 * 3)

# init <- list(center = truth$center + rnorm(d, 0, 1),
#              dispersion = 1,
#              Vt = truth$V0 + matrix(rnorm(q*d, 0, 1), nrow = d)
#              )
# init <- NULL

init_family <- function(x, weights, q, factor_family, sd_noise = 1){
      n = dim(x)[1]; d = dim(x)[2]
      mu <- family_initialize(x, weights, factor_family)
      eta <- factor_family$linkfun(mu)
      scale_eta <- scale(eta, scale = FALSE) # center
      center <- attr(scale_eta, "scaled:center") + rnorm(d, 0, sd_noise)
      S_ <- cov(eta)
      ss_ <- symm_eigen(S_)
      Vt <- sweep(ss_$vectors[, 1:q, drop = FALSE], 2, sqrt(ss_$values[1:q]), `*`) + matrix(rnorm(q*d, 0, sd_noise), nrow = d)
      dispersion = apply(weights * (x - mu)^2/factor_family$variance(mu), 2, mean)
return(list(center = center, dispersion = dispersion,Vt = Vt
))
}

init <- init_family(truth$X/truth$weights, truth$weights, q, factor_family, sd_noise = 0.5)



#load_dir = '/projectnb/dmfgrp/efm/OptiResult0108/'


if (algo_idx<=2){

  for (sample_idx in 1:3){
      sample_control <- list(sample_size = sample_list[sample_idx],
                             eval_size = 1500)


      save_name <- paste(load_dir, paste( algo_list[algo_idx],
                       name_list[family_idx], paste('s', sample_control$sample_size, sep =''),
                            paste('d', d, sep = ''),
                            paste('q', q, sep = ''),
                            paste('T', adam_control$max_epoch, sep= ''), sep = '_'),
                     '.RData', sep ='')


      #if (!file.exists(save_name)){
        start = Sys.time()
        efm_result <- efm(x = truth$X/factor_weights, lambda_prior = L_prior,
                          factor_family = factor_family,
                          rank = q, weights = factor_weights,
                          algo = algo_list[algo_idx],
                          start = init,
                          sample_control = sample_control,
                          adam_control = adam_control,
                          em_control = em_control,
                          eval_likeli = TRUE)
        end = Sys.time()
        efm_time = as.numeric(difftime(end, start, units = 's')) - efm_result$eval_time
        efm_result$efm_time = efm_time
        save(efm_result, file = save_name)
      #  }
  }# end of sample idx

}else{
  sample_idx = 1
  sample_control <- list(sample_size = sample_list[sample_idx],
                         eval_size = 1500)

  save_name <- paste(load_dir, paste( algo_list[algo_idx],
            name_list[family_idx], paste('s', sample_control$sample_size, sep =''),
            paste('d', d, sep = ''),
            paste('q', q, sep = ''),
            paste('T', adam_control$max_epoch, sep= ''), sep = '_'),
     '.RData', sep ='')
  #if (!file.exists(save_name)){
    start = Sys.time()
    efm_result <- efm(x = truth$X/factor_weights, lambda_prior = L_prior,
                      factor_family = factor_family,
                      rank = q, weights = factor_weights,
                      algo = algo_list[algo_idx],
                      start = init,
                      sample_control = sample_control,
                      adam_control = adam_control,
                      em_control = em_control,
                      eval_likeli = TRUE)
    end = Sys.time()
    efm_time = as.numeric(difftime(end, start, units = 's')) - efm_result$eval_time
    efm_result$efm_time = efm_time
    save(efm_result, file = save_name)
    #}
}

