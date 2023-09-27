setwd('/projectnb/dmfgrp/efm/')
source('./R/utils.R')
source('./R/efm.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

# [4.1 simulated data]
set.seed(1)
d = 10
n = 100
q = 2
dispersion_star = 2
adam_control = adam.control(max_epoch = 20, batch_size = 128,
                            step_size = 0.2, rho = 0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control = sample.control(sample_size = 500, eval_size = 500)

# factor_family = negative.binomial(dispersion_star)
factor_family = quasipoisson()
factor_weights = 1

# factor_family = Gamma(link = 'log')
# factor_weights = 1


# factor_family = quasibinomial()
# factor_family = binomial()
# factor_weights = matrix(rpois(n*d, 10) + 1, nrow =n )
#factor_weights = matrix(rep(10, n*d), nrow =n )


# [generate data]
prior_star = list(mean = rep(0, q),
     precision = rep(3, q))

V_star = matrix(rnorm(q*d, 0, 0.5), nrow = d)
#L_star = matrix(rnorm(n*q, 0, 0.5), nrow = n)
#L_star = matrix(rnorm(n*q, 0, 0.5), nrow = n)
L_star <- mvtnorm::rmvnorm(n, mean = prior_star$mean,
                             sigma = diag(1/prior_star$precision),
                             checkSymmetry = TRUE)


svd_star = svd(V_star, nu = q, nv = q)
negative_flag = c(1:q)[svd_star$u[1,]<0]
svd_star$u = efm_sign(svd_star$u)
svd_star$v = efm_sign(svd_star$v)
V_star = tcrossprod(svd_star$u, diag(svd_star$d))
L_star = tcrossprod(L_star, svd_star$v)

center_star = matrix(0, ncol = d)
eta_star = tcrossprod(L_star, V_star)
eta_star = sweep(eta_star, 2, center_star, '+')
X = generate_data(family= factor_family, eta = eta_star,
                  dispersion = dispersion_star, weights = factor_weights)
plot(density(X))


Vstart = svd(X, nu = q, nv = q)$v

center_start = center_star + rnorm(d, 0, 1)
dispersion_start = rnorm(d, 0, dispersion_star/3)
dispersion_start = dispersion_start - min(dispersion_start) + 0.1

# [initialize at the true]
# Vstart = V_star
# center_start = center_star
# dispersion_start = rep(dispersion_star, d)


true_likeli <- SML_neglikeli(V_star, factor_family, X, center = center_star, sample_control$eval_size,
                             dispersion = dispersion_star, weights = factor_weights, print_likeli = TRUE)


# [EM Newton]
# ngq <- 15
# control <- list(maxit = 20, epsilon = 1e-6, trace = TRUE)
# res <- fa_gqem(X, q, ngq, family =factor_family, control = control, Phi = dispersion_start)
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

init <- list(Vt = Vstart, center = center_star, dispersion = dispersion_start)
# init <- list(center = truth$center + rnorm(d, 0, 1),
#              dispersion = rnorm(d, 0, truth$phi/3),
#              Vt = svd(truth$X, nu = q, nv = q)$v)

# [Gradient]
lapl_result <- efm(X/factor_weights, factor_family,
                   rank = q, start = init,
                  factor_weights, algo = 'lapl',
                  adam_control = adam_control,
                  sample_control = sample_control, eval_likeli = TRUE,
                  lambda_prior = prior_star)



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




