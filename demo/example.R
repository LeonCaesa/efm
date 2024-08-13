library(efm)
#source('../R/utils.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

# [4.1 simulated data]
set.seed(2)

n = 500
q = 3
dispersion_star = 1
adam_control = adam.control(max_epoch = 20, batch_size = 128,
                            step_size = 0.2, rho = 0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control = sample.control(sample_size = 500, eval_size = 1500)

# [Poisson errors at the 5th step]
factor_family = poisson()
d = 10
factor_weights = 1


# [Binomial has its loss increasing]
# factor_family = binomial()
# d = 512
# factor_weights = matrix(rpois(n*d, 10) + 1, nrow =n )



# [generate data]
set.seed(2)
L_prior <- list(mean = rep(0, q),
                precision =  rep(1,q))

V_prior <- list(mean = rep(0, q),
                sigma = rep(1, q))
center_star <- rep(2, d)

truth <- generate_cov(n, d, L_prior, V_prior, center_star, family = factor_family,
                      weights = factor_weights, phi = dispersion_star)

plot(density(truth$X))




Vstart = svd(truth$X, nu = q, nv = q)$v
center_start = truth$center + rnorm(d, 0, 1)/10
dispersion_start = rnorm(d, 0, dispersion_star/3)
dispersion_start = dispersion_start - min(dispersion_start) + 0.1

# [initialize]
init_family <- function(x, weights, q, factor_family, sd_noise = 0) {
  n = dim(x)[1]; d = dim(x)[2]
  mu <- family_initialize(x, weights, factor_family)
  eta <- factor_family$linkfun(mu)
  scale_eta <- scale(eta, scale = FALSE) # center
  center <- attr(scale_eta, "scaled:center") + rnorm(d, 0, sd_noise)
  S_ <- cov(scale_eta)
  ss_ <- symm_eigen(S_)
  Vt <- sweep(ss_$vectors[, 1:q, drop = FALSE], 2, sqrt(ss_$values[1:q]), `*`) +
    matrix(rnorm(q*d, 0, sd_noise), nrow = d)
  dispersion = apply(weights * (x - mu) ^ 2 / factor_family$variance(mu), 2, mean)
  list(center = center, dispersion = dispersion, Vt = Vt)
}

#init <- init_family(truth$X/truth$weights, truth$weights, q, factor_family, sd_noise = 1)
init <- init_family(truth$X/truth$weights, truth$weights, q, factor_family, sd_noise = 0)




ngq <- 15
control <- list(maxit = 20, epsilon = 1e-6, trace = TRUE)

# [EM]
efm_result <- efm(truth$X/factor_weights, factor_family,
                                 rank = q, weights = factor_weights,
                                 start = init, algo ='em',  adam_control = adam_control,
                                 sample_control = sample_control, em_control = control,
                                 eval_likeli = TRUE, ngq= ngq,
                                 lambda_prior = L_prior)

