source('./R/utils.R')
source('./R/efm.R')
library(MASS)
library(matrixStats)
library(glmnet)
# [4.1 simulated data]
set.seed(1)
d = 5
n = 1000
q = 2
dispersion = 5
adam_control = adam.control(max_epoch = 30, batch_size = 128,
                            step_size = 0.5, rho = 0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control = sample.control(sample_size = 500, eval_size = 500)

#factor_family = negative.binomial(dispersion)
factor_family = quasipoisson()
factor_weights = 1


# factor_family = binomial()
# factor_weights = matrix(rpois(n*d, 10) + 1, nrow =n )
#factor_weights = matrix(rep(10, n*d), nrow =n )


# [generate data]
V_star = matrix(rnorm(q*d, 0, 1), nrow = d)
L_star = matrix(rnorm(n*q, 0, 1), nrow = n)
svd_star = svd(V_star, nu = q, nv = q)
negative_flag = c(1:q)[svd_star$u[1,]<0]
svd_star$u = efm_sign(svd_star$u)
svd_star$v = efm_sign(svd_star$v)
V_star = tcrossprod(svd_star$u, diag(svd_star$d))
L_star = tcrossprod(L_star, svd_star$v)

center_star = matrix(0, ncol = d)
eta_star = tcrossprod(L_star, V_star)
eta_star = sweep(eta_star, 2, center_star, '+')
X = generate_data(family= factor_family, eta = eta_star, dispersion = dispersion, weights = factor_weights)

Vstart = svd(X, nu = q, nv = q)$v;
center_start = center_star
dispersion_start = c(1,2,3,4,8)
# dispersion_start = rep(10, d)
lapl_result = efm(X/factor_weights, factor_family, center = center_start, rank = q, factor_weights, algo = 'lapl', start = V_star,
                  adam_control = adam_control, dispersion = dispersion_start,
                  sample_control = sample_control, eval_likeli = TRUE)
#lapl_result = efm_batch(Vstart, adam_control$batch_size, adam_control$step_size, X , factor_family, algo = 'lapl', eval_likeli = TRUE, eval_size =sample_control$sample_size)
print('done')




