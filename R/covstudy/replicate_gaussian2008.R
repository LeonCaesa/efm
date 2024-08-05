setwd('/projectnb/dmfgrp/efm/')
source('./R/efm.R')
source('./R/utils.R')
library(mvtnorm)
library(matrixStats)
library(MASS)
library(tidyverse)

metric_cov <- function(esti_cov, actual_cov){
  d <- dim(esti_cov)[1]
  eigen_frobenius <- eigen(esti_cov - actual_cov, only.values = TRUE)
  l2_frobenius <- sqrt(sum(eigen_frobenius$values^2))

  product = esti_cov %*% ginv(actual_cov)
  l1_entropy <- sum(diag(product)) - log (det(product)) - d

  sigma_half <- chol(ginv(actual_cov))
  inside_norm <- sigma_half %*% (esti_cov - actual_cov) %*% t(sigma_half)
  eigen_normalized <- eigen(inside_norm, only.values = TRUE)
  l2_normalized = sqrt(sum( eigen_normalized$values^2) ) *d^(-1/2)
  return(list(l2_frobenius, l1_entropy, l2_normalized))
}





q <- 3 ;n <- 756
L_prior <- list(mean = c(0.023558, 0.012989, 0.020714),
                precision =  1/ c(1.2507, 0.31564, 0.19303))

V_prior <- list(mean = c(0.78282, 0.51803, 0.41003),
                sigma = matrix( c(c(0.029145, 0.023873, 0.010184),
                                  c(0.023873, 0.053951, -0.006967),
                                  c(0.010184, -0.006967, 0.086856)), nrow = 3))
phi_prior<-list(alpah = 4.0713, beta= 0.1623)


factor_family= gaussian()
factor_weights = 1
# # [Optimization configuration]
# n_gq <- 15
# adam_control <- adam.control(max_epoch = 20, batch_size = 128,
#                              step_size = 0.2, rho = 0, abs_tol = 1e-8,
#                              beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8)
# sample_control <- sample.control(sample_size = 500, eval_size = 500)
# fagqem_control <- list(maxit = adam_control$max_epoch * as.integer(n / adam_control$batch_size),
#                        epsilon = 1e-8, trace = FALSE)

error_matrix <- matrix(ncol = 5)
metric_names = c('l2frobenius', 'l2entrophy', 'l2normalized')
d_list = seq(16, 1000, 100)

for (repeat_idx in 1:100){
  # fix factors
  for (d in d_list){
      print(c(repeat_idx, d))
      # L0 <- rmvnorm(n, mean = L_prior$mean,
      #             sigma = diag(1/L_prior$precision),
      #             checkSymmetry = TRUE)
      center_star <- rep(0, d)
      dispersion_star <- rgamma(d, shape = phi_prior$alpah, scale =phi_prior$beta)^2

      truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                            weights = factor_weights, phi = dispersion_star)#, L0 = L0)

      cov_true = truth$V0 %*% diag(1/L_prior$precision) %*% t(truth$V0) + diag(dispersion_star)
      cov_naive = cov(truth$X)

      c(l2fnorm, l2entrophy, l2normalized):=metric_cov(cov_naive, cov_true)
      output_metrics <- cbind(rbind(l2fnorm,l2entrophy,l2normalized), metric_names, d, 'naive', repeat_idx)
      error_matrix <- rbind(error_matrix, output_metrics)


      left_norm <- crossprod(truth$L0); right_norm <- crossprod(truth$L0, truth$X)
      V_hat <- t(backsolve(left_norm, right_norm))
      cov_common <- V_hat %*% cov(truth$L0) %*% t(V_hat)
      cov_single <- crossprod(truth$X - tcrossprod(truth$L0, V_hat))/n
      cov_factor08 <- cov_common + diag(diag(cov_single))
      c(l2fnorm, l2entrophy, l2normalized):=metric_cov(cov_factor08, cov_true)
      output_metrics <- cbind(rbind(l2fnorm,l2entrophy,l2normalized), metric_names, d, 'factor09', repeat_idx)
      error_matrix <- rbind(error_matrix, output_metrics)
  } # end of d loop
}# end of repeat loop


error_matrix <- error_matrix[-1, ]
error_df <- data.frame(error_matrix)
rownames(error_df ) <- 1:nrow(error_df )
colnames(error_df) =c('error', 'error_type', 'd', 'esti_method', 'repeat_idx')
error_df[c('error', 'd' )] = apply(error_df[c('error', 'd' )], 2, as.numeric)




plot_error = filter(error_df, esti_method %in% c('naive', 'factor09'))
ggplot(plot_error) + geom_boxplot(aes(x = as.factor(d), y = error,
                                      color = esti_method)) +
  facet_wrap(~error_type, scales = "free")


aggerror_df <- plot_error %>% group_by(d, esti_method, error_type) %>%
  summarise(avg_error = mean(error), sd_error = sd(error))

ggplot(aggerror_df) + geom_line(aes(x = d, y = avg_error, colour = esti_method, linetype = esti_method)) +
  facet_wrap(~error_type, scales = "free")

ggplot(aggerror_df) + geom_line(aes(x = d, y = sd_error, colour = esti_method, linetype = esti_method)) +
  facet_wrap(~error_type, scales = "free")
