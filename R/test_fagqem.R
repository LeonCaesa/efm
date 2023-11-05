setwd('/Users/caesa/Desktop/efm/')
source('./R/efm.R')
source('./R/utils.R')
library(mvtnorm)
library(matrixStats)
library(MASS)

# [4.1 simulated data]
set.seed(1)
d = 10
n = 50
q = 3


# dispersion_star = 20
# factor_family = negative.binomial(dispersion_star)
# factor_weights = 1

# dispersion_star = 2
# factor_family = quasipoisson()
# factor_weights = 1

dispersion_star = 1
factor_family = binomial()
factor_weights = 10

# [Generate data]
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
init <- list(Vt= Vstart, center = center_start, dispersion = dispersion_start)



# [Speed and Binomial test]
ngq <- 15
control <- list(maxit = 20, epsilon = 1e-6, trace = TRUE)

# [Luis's code with dispersion support -- Faster]
start = Sys.time()
res2 <- fa_gqem2(truth$X/factor_weights, q, ngq, family = factor_family, control = control,
                 eval_likeli = TRUE, weights = factor_weights, Phi = dispersion_start)
end = Sys.time()
time2 = end - start
print(time2)

# [Original fagqem code --Slower]
start = Sys.time()
res <- fa_gqem(truth$X/factor_weights, q, ngq, family =factor_family, control = control,
               eval_likeli = TRUE, weights = factor_weights, Phi = dispersion_start)
end = Sys.time()
time1 = end - start
print(time1)
