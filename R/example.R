setwd('/projectnb/dmfgrp/efm/')
source('./R/utils.R')
source('./R/efm.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

# [4.1 simulated data]
set.seed(1)
d = 50
n = 500
q = 3
dispersion_star = 20
adam_control = adam.control(max_epoch = 20, batch_size = 128,
                            step_size = 0.2, rho = 0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control = sample.control(sample_size = 500, eval_size = 500)

factor_family = negative.binomial(dispersion_star)
# factor_family = quasipoisson()
factor_weights = 1

# factor_family = Gamma(link = 'log')
# factor_weights = 1


# factor_family = quasibinomial()
# factor_family = binomial()
# factor_weights = matrix(rpois(n*d, 10) + 1, nrow =n )
#factor_weights = matrix(rep(10, n*d), nrow =n )


# [generate data]
L_prior <- list(mean = c(0, 0, 0),
                precision =  1/ c(1, 1, 1))

V_prior <- list(mean = c(0, 0 ,0),
                sigma = c(0.5, 0.5, 0.5))
center_star <- center_star <- rep(0, d)


truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                      weights = factor_weights, phi = dispersion_star)

plot(density(truth$X))




Vstart = svd(truth$X, nu = q, nv = q)$v
center_start = truth$center + rnorm(d, 0, 1)/10
dispersion_start = rnorm(d, 0, dispersion_star/3)
dispersion_start = dispersion_start - min(dispersion_start) + 0.1

# [initialize at the true]
# Vstart = V_star
# center_start = center_star
# dispersion_start = rep(dispersion_star, d)

init <- NULL
# init <- list(Vt = Vstart, center = center_start, dispersion = dispersion_start)
# init <- list(center = truth$center + rnorm(d, 0, 1),
#              dispersion = rnorm(d, 0, truth$phi/3),
#              Vt = svd(truth$X, nu = q, nv = q)$v)



# [Gradient]
lapl_result <- efm(truth$X/factor_weights, factor_family,
                                 rank = q, weights = factor_weights,
                                 start = init, algo ='lapl',  adam_control = adam_control,
                                 sample_control = sample_control, eval_likeli = TRUE,
                                 lambda_prior = L_prior)


# [EM Newton]
ngq <- 15
control <- list(maxit = 20, epsilon = 1e-6, trace = TRUE)
res <- fa_gqem(truth$X, q, ngq, family = factor_family, lambda_prior= L_prior,
               start = init, control = control, eval_likeli = TRUE,
               eval_size = sample_control$eval_size)



#
# cbind(res$alpha, res$V, c(center_star), V_star)
# cov(X); cor(X)
# (efm_esti_cov <- marg_var(res$alpha, res$V, res$family, ngq)); cov2cor(cx)
# efm_true_cov <- marg_var(center_star, V_star, res$family, ngq)
#
# mse_efm = mean((efm_esti_cov- efm_true_cov)^2)
# mse_num = mean((cov(X)- efm_true_cov)^2)
# print(mse_efm)
# print(mse_num)



# [Param comparison]
plot(lapl_result$like_list)
cbind(dispersion_start, lapl_result$dispersion, dispersion_star)
matrix(cbind(center_start, lapl_result$center, center_star), ncol = 3)
cbind(Vstart, lapl_result$V, V_star)


# [Covariance comparison]
num_esti_cov = cov(X)
efm_esti_cov = marg_var(lapl_result$center, lapl_result$V, lapl_result$family, ngq = 15)
efm_true_cov = marg_var(center_star, V_star, lapl_result$family, ngq = 15)


mse_efm = mean((efm_esti_cov- efm_true_cov)^2)
mse_num = mean((num_esti_cov- efm_true_cov)^2)
print(mse_efm)
print(mse_num)




