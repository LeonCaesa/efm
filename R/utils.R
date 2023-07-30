norm2 <- function (x) norm(as.matrix(x[!is.na(x)]), "2")
normalize <- function (x, margin = 2)
  sweep(x, margin, apply(x, margin, norm2), `/`)

pdf_calc <-
  function(family,
           weights,
           dispersion = 1,
           log_ = FALSE) {
    if (grepl('^Negative Binomial', family$family))
      (family$family = 'negbinom')
    family_pdf <- switch(
      family$family,
      gaussian = function(x, mu, weights)
        weights * dnorm(x, mu,  dispersion, log = log_),
      poisson = function(x, mu, weights)
        weights *  dpois(x, mu, log = log_),
      binomial = function(x, mu, weights)
        weights * dbinom(x * weights, weights, mu, log = log_),
      negbinom = function(x, mu, weights)
        weights *  dnbinom(x,  dispersion, mu = mu, log = log_),
      Gamma = function(x, mu, weights)
        weights *  dgamma(
          x,
          shape =  dispersion,
          scale = mu / dispersion,
          log = log_
        ),
      stop("Family `", family$family, "` not recognized")
    )
    return (family_pdf)
  }


gsym_solve <- function (A, b, tol = sqrt(.Machine$double.eps)) {
  ea <- eigen(A, symmetric = TRUE)
  V <- ea$vectors; d <- ea$values
  valid <- d > max(tol * d[1], 0)
  if (!all(valid)) {
    V <- V[, valid, drop = FALSE]; d <- d[valid]
  }
  V %*% sweep(crossprod(V, b), 1, d, `/`)
}



comput_mupos = function(L_row, Vt, factor_family, scale_weights = 1) {
  q = dim(Vt)[2]
  d = dim(Vt)[1]
  Weight_row = tcrossprod(L_row, Vt)
  mu_hat = c(factor_family$linkinv(Weight_row))

  diag_term = diag(scale_weights * factor_family$variance(mu_hat), nrow = d)
  hessian = tcrossprod(crossprod(Vt, diag_term), t(Vt))
  mu_pos = crossprod(ginv(ginv(hessian) + diag(1, nrow= q)), L_row)
  return( mu_pos)
}


comput_CholSigma = function(L_row, Vt, factor_family, scale_weights = 1){
  q = dim(Vt)[2]
  d = dim(Vt)[1]
  mu_hat = c(factor_family$linkinv(tcrossprod(L_row, Vt)))
  diag_term = diag(scale_weights * factor_family$variance(mu_hat), nrow = d)
  hessian = tcrossprod(crossprod(Vt, diag_term), t(Vt))
  Sigma_pos = ginv(diag(1, nrow= q) + hessian) + diag(1e-08, nrow = q)
  Sigma_chol = chol(Sigma_pos)
  return(Sigma_chol)
}


simu_pos = function(mu_pos, CholSigma, sample_size = 1){
  q = length(mu_pos)
  simu_temp = tcrossprod(CholSigma, matrix(rnorm(sample_size * q), nrow = sample_size))
  return( sweep(simu_temp, 1, mu_pos, '+') )
}


# [gradient calculation --- Laplacian]
ridge_coef <- function(X_vec, weight_vec, Vt, center, factor_family){
  d = dim(Vt)[1]
  sd_scalar = sqrt(var(X_vec)*(d-1)/d)
  pen_result <- glmnet(x = Vt, y= X_vec, family = factor_family, alpha = 0, lambda = 1,
                       weights = weight_vec, offset = center,
                       intercept = FALSE, standardize = FALSE, thresh= 1e-10,
                       type.logistic = c("Newton"))
  #print(Vt)
  return(as.vector(coef(pen_result, s = sd_scalar * 1/d, exact = TRUE, x = Vt, y = X_vec,
                 family = factor_family, offset = center,
                 weights = weight_vec))[-1])
}



Lapl_grad <- function(X_batch, Vt, factor_family, weights, dispersion, center, q = dim(Vt)[2]){
  n <- dim(X_batch)[1]
  d <- dim(X_batch)[2]
  Vt <- matrix(Vt, nrow = d, ncol =q)

  L_mle = t(mapply(ridge_coef, X_vec = asplit(X_batch, 1),
                   weight_vec = asplit(weights, 1),
                   MoreArgs = list(Vt = Vt, factor_family = factor_family, center = c(center))))

  eta_mle = sweep(tcrossprod(L_mle, Vt), 2 , center, '+')
  mu_mle = factor_family$linkinv(eta_mle)
  g = factor_family$mu.eta(eta_mle); dim(g) = dim(eta_mle)
  v = factor_family$variance(mu_mle); dim(v) = dim(eta_mle)
  scales = (mu_mle - X_batch) * weights * g/v # to do: confirm the weights
  grad_V = 1/ dispersion * t(crossprod(scales, L_mle))

  grad_center = colSums(scales/dispersion)


  return(list(grad_V = grad_V, grad_center = grad_center))
}

# [gradient calculation --- SML]

SML_grad <- function(Vt, factor_family, X, q, sample_size, dispersion = 1, weights = 1){
  n <- dim(X)[1]
  d <- dim(X)[2]
  Vt <- matrix(Vt, nrow= d, ncol =q)
  family_pdf = pdf_calc(factor_family, weights = weights, dispersion = dispersion, log_ = TRUE)

  likeli_simu <- array(dim = c(n, d, sample_size))
  grad_simu <- array(dim = c(n, d, q, sample_size))
  for (b in 1:sample_size){
    L_sample <-  matrix(rnorm(n * q), nrow = n)
    eta_simu <- tcrossprod(L_sample, Vt)
    mu_simu <- factor_family$linkinv(eta_simu)
    variance_simu <- factor_family$variance(mu_simu)
    mueta_simu <- factor_family$mu.eta(eta_simu)
    grad_scale = 1/dispersion * (mu_simu- X) * weights * mueta_simu/variance_simu# to do: check glm_weight multiply
    for( j in 1:d){
      for(l in 1:q){
        grad_simu[,j, l, b]= grad_scale[,j] * L_sample[,l]}  # end of q iter
    }  # end of d iter
    likeli_simu[,,b] <- family_pdf(X, mu = mu_simu, weights = weights)
  } # end of sample iter
  mc_likeli = exp(apply(likeli_simu, c(1, 2), logSumExp))/sample_size
  mc_likeli = replicate(q, mc_likeli)
  mc_grad = apply(grad_simu, c(1,2,3), mean)

  total_grad = apply(mc_grad/mc_likeli, c(2,3), sum)
  #to avoid numerical issue of dividing a small probability
  total_grad[is.na(total_grad)] = 0 #todo: discuss this numerical workaround with Luis
  return( t(total_grad))
}



# [gradient calculation --- Posterior Sampling]

PosSample_Moments <-function(LX_row, Vt, factor_family, dispersion, weights){
  d = dim(Vt)[1]
  q = dim(Vt)[2]
  L_row = matrix(LX_row[1:q], nrow = 1)
  X_row = matrix(LX_row[-c(1:q)], nrow = 1)
  eta_row = tcrossprod(L_row, Vt)
  mu_row =  factor_family$linkinv(eta_row)
  var_row = factor_family$variance(mu_row)

  var_row = sweep(var_row, 2, weights , '*')
  hessian_L = 1/dispersion * crossprod(sweep(Vt, 1, var_row, "*"), Vt)

  L_pos = matrix(tcrossprod(ginv(ginv(hessian_L) + diag(1, nrow= q)), L_row), ncol=1)
  Sigma_pos = ginv(diag(1, nrow= q) + hessian_L) + diag(1e-08, nrow = q)

  return(list(L_pos = L_pos, Sigma_pos = Sigma_pos))
}

PosSample_GradRow <- function(L_row, X_row, weight_row, Vt, factor_family, dispersion, sample_size){

  d = dim(Vt)[1];q = dim(Vt)[2]
  LX_row = cbind(L_row, X_row)
  Pos_moments <- PosSample_Moments(LX_row, Vt, factor_family, dispersion, weight_row)
  L_sample <- PosSample_Draw(Pos_moments, sample_size)

  eta_pos = tcrossprod(L_sample, Vt)
  mu_pos = factor_family$linkinv(eta_pos)
  var_pos = factor_family$variance(mu_pos)
  mueta_pos = factor_family$mu.eta(eta_pos)


  scales = sweep(mu_pos, 2, X_row, '-')
  scales = scales * mueta_pos/var_pos
  scales = sweep(scales, 2, weight_row, '*') # [todo: examine glm weights]
  grad_return = 1/dispersion * crossprod(L_sample, scales)/sample_size
  return( c(grad_return))
}

PosSample_Draw <- function(Pos_moments, sample_size){
  L_pos = Pos_moments$L_pos
  Sigma_chol = chol(Pos_moments$Sigma_pos)

  sd_normal = matrix(rnorm(sample_size * dim(L_pos)[1]), nrow = sample_size)
  L_sample = t(tcrossprod(Sigma_chol, sd_normal))
  L_sample = sweep(L_sample, 2, L_pos,'+')
  return(L_sample)
}


PosSample_grad <-function(X_batch, Vt, factor_family, dispersion, weights,
                          sample_size = sample_size, q= dim(Vt)[2]){
  n <- dim(X_batch)[1];d <- dim(X_batch)[2];q = dim(Vt)[2]
  Vt <- matrix(Vt, nrow= d, ncol =q)

  L_mle = batch_mle(X_batch, Vt, factor_family, q, weights)
  grad = mapply(PosSample_GradRow,
                L_row = asplit(L_mle, 1),
                X_row = asplit(X_batch, 1),
                weight_row = asplit(weights, 1),
                MoreArgs = list(Vt = Vt, factor_family = factor_family, dispersion = dispersion, sample_size = sample_size))
  return (matrix(rowSums(grad), nrow = q))
}




# [gradient calculation --- Fast Gaussian Plug In]
FastGaussian_grad <- function(X_batch, Vt, factor_family, dispersion, weights,  q= dim(Vt)[2]){
  n <- dim(X_batch)[1]
  d <- dim(X_batch)[2]
  Vt <- matrix(Vt, nrow= d, ncol =q)

  q = dim(Vt)[2]
  d = dim(Vt)[1]
  L_mle = batch_mle(X_batch, Vt, factor_family, q)
  LX_matrix = cbind(L_mle, X_batch)
  grad = apply(LX_matrix, 1, FastGaussian_GradRow, Vt =Vt,
               factor_family = factor_family, dispersion = dispersion,
               weights = weights)
  grad = matrix(rowSums(grad), ncol =d)
}



FastGaussian_GradRow <- function(LX_row, Vt, factor_family, dispersion, weights){


  q = dim(Vt)[2]
  # family_pdf = pdf_calc(factor_family, weights = weights, dispersion = dispersion, log_ = TRUE)

  L_row = matrix(LX_row[1:q], nrow = 1)
  X_row = matrix(LX_row[-c(1:q)], nrow = 1)
  eta_row = tcrossprod(L_row, Vt)
  mu_row =  factor_family$linkinv(eta_row)
  var_row = factor_family$variance(mu_row)
  mueta_row = factor_family$mu.eta(eta_row)
  hessian_L = 1/dispersion * crossprod(sweep(Vt, 1, var_row, "*"), Vt)

  L_pos = matrix(tcrossprod(ginv(ginv(hessian_L) + diag(1, nrow= q)), L_row), nrow=1)
  Sigma_pos = ginv(diag(1, nrow= q) + hessian_L) + diag(1e-08, nrow = q)

  eta_pos = tcrossprod(L_pos, Vt)
  mu_pos = factor_family$linkinv(eta_pos)
  # [todo: fix glm weights]
  scales = (mu_pos - X_row) * mueta_row/var_row
  # scales = (X_row - mu_pos) * mueta_row/var_row
  b_hat = 1/dispersion * crossprod(L_pos,  scales)

  # pos_likeli = sum(family_pdf(X_row, mu_pos, weights = weights)) + dmvnorm(L_pos, mean = rep(0, q), sigma = Sigma_pos, log = TRUE)
  #pos_likeli = dmvnorm(X_row, mean = mu_row, sigma =diag(c(var_row), nrow = d)) * dmvnorm(L_pos, mean = rep(0,q), sigma = Sigma_pos)
  # pos_likeli = dmvnorm(L_pos, mean = mu_pos, sigma = Sigma_pos)
  constant = (2 * pi)^(q/2) * det(Sigma_pos)^(1/2)
  return(b_hat * constant)
}
