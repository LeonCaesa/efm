norm2 <- function (x) norm(as.matrix(x[!is.na(x)]), "2")
normalize <- function (x, margin = 2)
  sweep(x, margin, apply(x, margin, norm2), `/`)

mat_mult <- function (A, b, mult_cond = is.vector(A))
  if (mult_cond) A * b else A %*% b

mat_add <- function (A, B, add_cond = is.vector(A)) {
  if (add_cond) diag(B) <- diag(B) + A else B <- B + A
  B
}

# "Robust" eigendecomposition of symmetric matrix A
symm_eigen <- function (A, rank_tol = sqrt(.Machine$double.eps)) {
  ea <- eigen(A, symmetric = TRUE)
  V <- ea$vectors; d <- ea$values
  valid <- ea$values > max(rank_tol * ea$values[1], 0)
  if (!all(valid)) {
    ea$vectors <- ea$vectors[, valid, drop = FALSE]
    ea$values <- ea$values[valid]
  }
  ea
}

# If a symmetric matrix A has eigen decomposition
# A = V * D * V' = (D_2 * V')' * (D_2 * V') = C' * C,
# where D_2 = D^(1/2) and C = D_2 * V',
# performs one of the following operations:
# solve = FALSE, transpose = FALSE: C * b = D_2 * V' * b
# solve = FALSE, transpose = TRUE: C' * b =  V * D_2 * b
# solve = TRUE, transpose = FALSE: C^(-1) * b = V * D_2^(-1) * b
# solve = TRUE, transpose = TRUE: C'^(-1) * b = D_2^(-1) * V' * b
symm_mult <- function (ea, b, solve = FALSE, transpose = FALSE) {
  d_2 <- sqrt(ea$values)
  if (!xor(solve, transpose)) # apply V' first?
    sweep(crossprod(ea$vectors, b), 1, d_2, ifelse(solve, `/`, `*`))
  else
    ea$vectors %*% sweep(b, 1, d_2, ifelse(solve, `/`, `*`))
}




#' @export
#' Var(X) = E_L[Var(X|L) + Var(E(X|L))]
marg_var <- function (center, V, family, ngq = 15, ...) {
  if (!is.vector(center)) center = c(center)
  gq <- gaussquadr::GaussQuad$new("gaussian", ngq, ncol(V), ...)
  fam <- if (is.function(family)) family() else family
  ev <- gq$E(function (l) fam$variance(fam$linkinv(center + V %*% l)))
  mv <- gq$E(function (l) fam$linkinv(center + V %*% l))
  ve <- gq$E(function (l) tcrossprod(fam$linkinv(center + V %*% l) - mv))
  matrix(ve, nrow = nrow(V)) + diag(ev) # question: why add diagonal?
}


check_canonical <- function (fam) {
  ((!is.null(fam$linkcan)) && fam$linkcan) ||
    ((fam$family == "poisson" || fam$family == "quasipoisson") &&
       fam$link == "log") ||
    ((fam$family == "binomial" || fam$family == "quasibinomial") &&
       fam$link == "logit") ||
    (fam$family == "gaussian" && fam$link == "identity")
}


dquasipois <- function(x, mu, dispersion, log = TRUE){
  var <- dispersion * mu
  #add_const <- x - x * log(x)
  prob_log = 1/dispersion * (x * log(mu) - mu) - 0.5 * log(2 * pi * var) #+ add_const
  if (log){return(prob_log)
  }else{return(exp(prob_log))}
}

dquasibinom <- function(x, mu, dispersion, log = TRUE){
  var <- dispersion * mu * (1-mu)
  logits <- x * log(mu/(1-mu)) + log(1-mu)
  prob_log = 1/dispersion * logits - 0.5 * log(2 * pi * var) #+ add_const
  if (log){return(prob_log)
  }else{return(exp(prob_log))}
}


pdf_calc <-
  function(family,
           weights,
           dispersion = 1,
           log_ = FALSE) {
    if (grepl('^Negative Binomial', family$family))
      (family$family = 'negbinom')
    family_pdf <- switch(
      family$family,
      gaussian = function(x, mu, weights, dispersion)
        weights * dnorm(x, mu,  dispersion, log = log_),
      poisson = function(x, mu, weights, dispersion = 1)
        weights *  dpois(x, mu, log = log_),
      quasipoisson = function(x, mu, weights, dispersion)
        weights * dquasipois(x, mu, dispersion, log = log_),
      quasibinomial= function(x, mu, weights, dispersion)
        weights * dquasibinom(x, mu, dispersion, log = log_),
      binomial = function(x, mu, weights, dispersion = 1)
        weights * dbinom(x * weights, weights, mu, log = log_),
      negbinom = function(x, mu, weights, dispersion = 1)
        weights *  dnbinom(x,  1/(factor_family$variance(1)-1), mu = mu, log = log_),
      Gamma = function(x, mu, weights, dispersion)
        weights *  dgamma(
          x,
          shape =  1/dispersion,
          scale = mu * dispersion,
          log = log_
        ),
      stop("Family `", family$family, "` not recognized")
    )
    return (family_pdf)
  }

check_DispersionUpdate <- function(family){
    update_families = c('quasipoisson', 'quasibinomial', 'Gaussian', 'Gamma')
    if (family$family %in% update_families){
      return (TRUE)
    }
    return (FALSE)
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


# [method for experiment generation]
# todo: add prior support
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
  simu_temp = tcrossprod(CholSigma, matrix(rnorm(sample_size * q),
                                           nrow = sample_size))
  return( sweep(simu_temp, 1, mu_pos, '+') )
}



bsglm <- function (x, y, prior_coef, weights = NULL, offset = NULL,
                   start = NULL, family = gaussian(), dispersion = 1,
                   control = list(), return_hessian = FALSE) {
  control <- do.call("glm.control", control)
  nobs <- nrow(x); nvars <- ncol(x)
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  # [ initialize ]
  if (!is.null(start)) {
    eta <- drop(if (nvars == 1L) x * start else x %*% start)
  } else {
    if (family$family == "gaussian") etastart <- 1
    eval(family$initialize)
    eta <- family$linkfun(mustart)
  }
  eta <- eta + offset; mu <- family$linkinv(eta)
  beta0 <- mat_mult(prior_coef$precision, drop(prior_coef$mean))
  dev <- Inf
  # [ iterate ]
  for (it in 1:control$maxit) {
    varmu <- family$variance(mu)
    W <- weights * varmu / dispersion
    residuals <- weights * (y - mu) / dispersion
    if (!check_canonical(family)) { # adjust weights?
      thetaeta <- family$mu.eta(eta) / varmu
      W <- W * thetaeta ^ 2
      residuals <- residuals * thetaeta
    }
    z <- crossprod(x, W * (eta - offset) + residuals) + beta0
    H <- mat_add(prior_coef$precision, crossprod(x, x * c(W))) # Hessian
    eh <- symm_eigen(H)
    #beta <- gsym_solve_bsglm(eh, z)
    beta <- gsym_solve(H, z)

    eta <- drop(mat_mult(x, beta, nvars == 1)) + offset
    mu <- family$linkinv(eta)
    bd <- beta - prior_coef$mean
    dev_new <- sum(family$dev.resids(y, mu, weights)/dispersion) +
      sum(bd * mat_mult(prior_coef$precision, bd))  # FIXME: check
    if (control$trace) message("<", it, "> dev = ", dev_new)
    if (it > 1 && abs((dev_new - dev) / (dev + .1)) < control$epsilon) break
    dev <- dev_new
  }
  mueta_square <- family$mu.eta(eta)^2
  varmu <- family$variance(mu) * dispersion
  diag_s = weights / varmu * mueta_square
  # S = n times n , non-diagonal matrix
  if (return_hessian){
    return (list(coef = beta, hessian = eh, deviance = dev, diag_s = diag_s))
  }else{
    return (beta)
  }
}

Lapl_grad <- function(X_batch, Vt,
                      factor_family, weights,
                      dispersion, center, q = dim(Vt)[2],
                      lambda_prior = list(mean = rep(0, q),
                                    precision = rep(1, q))
                      ){

  n <- dim(X_batch)[1]; d <- dim(X_batch)[2]
  Vt <- matrix(Vt, nrow = d, ncol =q)


  L_mle = t(mapply(bsglm, y = asplit(X_batch, 1),
           weights = asplit(weights, 1),
           MoreArgs = list(x = Vt, prior_coef = lambda_prior, dispersion = dispersion,
                           family = factor_family, offset = c(center))
           ))
  # question, returning a list of variable

  # end = Sys.time()
  # print(end- start)

  #browser()
  #head(L_mle)


  eta_mle <- sweep(tcrossprod(L_mle, Vt), 2 , center, '+')
  mu_mle <- factor_family$linkinv(eta_mle)
  g <- factor_family$mu.eta(eta_mle); dim(g) <- dim(eta_mle)
  v <- factor_family$variance(mu_mle); dim(v) <- dim(eta_mle)


  scales <- (mu_mle - X_batch) * weights * g/v
  scale_dispersion <-  1 - sweep(weights * (mu_mle - X_batch)^2 / v, 2, dispersion, '/')


  scales <- sweep(scales, 2, dispersion, '/')
  grad_V <- crossprod(scales, L_mle)
  grad_center <- colSums(scales)
  grad_dispersion <- colSums(sweep(scale_dispersion, 2, dispersion * 2, '/'))


  return(list(grad_V = grad_V, grad_center = grad_center,
              grad_dispersion = grad_dispersion))
}

# [gradient calculation --- SML]
SML_grad <- function(Vt, factor_family, X, q, center,
                     sample_size, dispersion = 1, weights = 1,
                     lambda_prior = list(mean = rep(0, q),
                                         precision = rep(1, q))
                     ){
  n <- dim(X)[1]; d <- dim(X)[2]
  Vt <- matrix(Vt, nrow = d, ncol = q)
  family_pdf <- pdf_calc(factor_family, weights = weights,
                         dispersion = dispersion, log_ = TRUE)

  likeli_simu <- array(dim = c(n, d, sample_size))
  grad_simu <- array(dim = c(n, d, q, sample_size))
  grad_simu_center <- array(dim = c(n, d, sample_size))
  grad_simu_dispersion <- matrix(nrow = sample_size, ncol = d)
  for (b in 1:sample_size){
    L_sample <- mvtnorm::rmvnorm(n, mean = lambda_prior$mean,
                                 sigma = diag(1/lambda_prior$precision),
                                 checkSymmetry = TRUE)
    eta_simu <- tcrossprod(L_sample, Vt)
    eta_simu <- sweep(eta_simu, 2, center, '+')
    mu_simu <- factor_family$linkinv(eta_simu)
    variance_simu <- factor_family$variance(mu_simu)
    mueta_simu <- factor_family$mu.eta(eta_simu)

    grad_scale <- (mu_simu- X) * weights * mueta_simu/variance_simu
    grad_scale <- sweep(grad_scale, 2, dispersion, '/' )

    grad_Wpearson <- (mu_simu- X)^2 /variance_simu * weights
    grad_disperscale <- 1 - sweep(grad_Wpearson, 2, dispersion, '/')
    grad_simu_dispersion[b,] <- colSums(sweep(grad_disperscale, 2, dispersion * 2, '/'))


    for( j in 1:d){
      for(l in 1:q){
        grad_simu[,j, l, b] <- grad_scale[,j] * L_sample[,l]}  # end of q iter
    }  # end of d iter
    grad_simu_center[,,b] <- grad_scale
    likeli_simu[,,b] <- family_pdf(X, mu = mu_simu, dispersion = dispersion, weights = weights)
  } # end of sample iter
  mc_likeli <- exp(apply(likeli_simu, c(1, 2), logSumExp))/sample_size
  mc_likeli <- replicate(q, mc_likeli)
  mc_grad_v <- apply(grad_simu, c(1,2,3), mean)
  mc_grad_center <- apply(grad_simu_center, c(1,2), mean)

  total_grad_v <- apply(mc_grad_v/mc_likeli, c(2,3), sum)
  total_grad_center <- colSums(mc_grad_center/mc_likeli[,,1])

  #to avoid numerical issue of dividing a small probability
  total_grad_v[is.na(total_grad_v)] <- 0 #todo: discuss this numerical workaround with Luis
  total_grad_center[is.na(total_grad_center)] <- 0

  total_grad_dispersion <- colMeans(grad_simu_dispersion)
  return( list(grad_V = total_grad_v,
               grad_center = total_grad_center,
               grad_dispersion = total_grad_dispersion
               ))
}



# [gradient calculation --- Posterior Sampling]
PosSample_Moments <-function(LX_row, Vt, center, factor_family,
                             dispersion, weights,
                             lambda_prior = list(mean = rep(0, q),
                                                 precision = rep(1, q))
                             ){
  d = dim(Vt)[1]; q = dim(Vt)[2]
  L_row = matrix(LX_row[1:q], nrow = 1)
  X_row = matrix(LX_row[-c(1:q)], nrow = 1)
  eta_row = tcrossprod(L_row, Vt)
  eta_row = sweep(eta_row, 2, center, '+')
  mu_row =  factor_family$linkinv(eta_row)
  var_row = factor_family$variance(mu_row)
  mueta_row = factor_family$mu.eta(eta_row)


  scale_var = 1/sweep(var_row, 2, dispersion/weights , '*')
  scale_var = sweep(scale_var, 2, mueta_row, '*')
  hessian_L = crossprod(sweep(Vt, 1, scale_var, "*"), Vt)
  inv_hessian_L = ginv(hessian_L)

  lambda_prior$Sigma = diag(1/lambda_prior$precision)
  hessian_add =  inv_hessian_L + lambda_prior$Sigma
  hessian_addinv = ginv(hessian_add)

  L_pos1 <- lambda_prior$Sigma %*% tcrossprod(hessian_addinv, L_row)
  L_pos2 <- inv_hessian_L %*% hessian_addinv %*% lambda_prior$mean
  L_pos = L_pos1 + L_pos2

  Sigma_pos = inv_hessian_L %*% hessian_addinv %*% lambda_prior$Sigma
  return(list(L_pos = L_pos, Sigma_pos = Sigma_pos))
}

# todo: PosSample_GradRowV, PosSample_GradRowPhi, PosSample_GradRowCenter
PosSample_GradRowV <- function(L_row, X_row, center,
                               weight_row, Vt, factor_family,
                               dispersion, sample_size,
                               lambda_prior = list(mean = rep(0, q),
                                                   precision = rep(1, q))
                               ){

  d = dim(Vt)[1];q = dim(Vt)[2]
  LX_row = c(L_row, X_row)

  Pos_moments <- PosSample_Moments(LX_row = LX_row, Vt = Vt,
                                   center = center,
                                   factor_family = factor_family,
                                   dispersion = dispersion,
                                   weights = weight_row,
                                   lambda_prior = lambda_prior)
  L_sample <- PosSample_Draw(Pos_moments, sample_size)

  eta_pos = tcrossprod(L_sample, Vt)
  eta_pos = sweep(eta_pos, 2, center, "+")
  mu_pos = factor_family$linkinv(eta_pos)
  var_pos = factor_family$variance(mu_pos)
  mueta_pos = factor_family$mu.eta(eta_pos)


  scales = sweep(mu_pos, 2, X_row, '-')
  scales = scales * mueta_pos/var_pos
  scales = sweep(scales, 2, weight_row/dispersion, '*')
  grad_v = crossprod(scales, L_sample)/sample_size
  return( grad_v)
}

PosSample_GradRowPhi <- function(L_row, X_row, center,
                                 weight_row, Vt, factor_family,
                                 dispersion, sample_size,
                                 lambda_prior = list(mean = rep(0, q),
                                                     precision = rep(1, q))
                                 ){

  d = dim(Vt)[1]; q = dim(Vt)[2]
  LX_row = c(L_row, X_row)
  Pos_moments <- PosSample_Moments(LX_row = LX_row, Vt = Vt,
                                   center = center,
                                   factor_family = factor_family,
                                   dispersion = dispersion,
                                   weights = weight_row,
                                   lambda_prior = lambda_prior)
  L_sample <- PosSample_Draw(Pos_moments, sample_size)

  eta_pos = tcrossprod(L_sample, Vt)
  eta_pos = sweep(eta_pos, 2, center, "+")
  mu_pos = factor_family$linkinv(eta_pos)
  var_pos = factor_family$variance(mu_pos)

  scales = -sweep(mu_pos, 2, X_row, '-')
  scales = 1 - sweep(scales^2, 2, dispersion/weight_row, '/')
  scales = sweep(scales, 2, dispersion, '/')
  grad_phi = colMeans(scales)
  return(c(grad_phi))
}

PosSample_GradRowCenter <- function(L_row, X_row, center,
                                    weight_row, Vt, factor_family,
                                    dispersion, sample_size,
                                    lambda_prior = list(mean = rep(0, q),
                                                        precision = rep(1, q))){

  d = dim(Vt)[1];q = dim(Vt)[2]
  LX_row = c(L_row, X_row)
  Pos_moments <- PosSample_Moments(LX_row = LX_row, Vt = Vt,
                                   center = center,
                                   factor_family = factor_family,
                                   dispersion = dispersion,
                                   weights = weight_row,
                                   lambda_prior = lambda_prior)
  L_sample <- PosSample_Draw(Pos_moments, sample_size)

  eta_pos = tcrossprod(L_sample, Vt)
  eta_pos = sweep(eta_pos, 2, center, "+")
  mu_pos = factor_family$linkinv(eta_pos)
  var_pos = factor_family$variance(mu_pos)
  mueta_pos = factor_family$mu.eta(eta_pos)
  scales = sweep(mu_pos, 2, X_row, '-')
  scales = scales * mueta_pos/var_pos
  scales = sweep(scales, 2, weight_row, '*')
  grad_center = 1/dispersion * colMeans(scales)
  return( grad_center)
}
PosSample_Draw <- function(Pos_moments, sample_size){
  L_pos = Pos_moments$L_pos
  Sigma_chol = chol(Pos_moments$Sigma_pos)

  sd_normal = matrix(rnorm(sample_size * dim(L_pos)[1]), nrow = sample_size)
  L_sample = t(tcrossprod(Sigma_chol, sd_normal))
  L_sample = sweep(L_sample, 2, L_pos,'+')
  return(L_sample)
}


PosSample_grad <-function(X_batch, Vt, center, factor_family,
                          dispersion, weights, sample_size = sample_size,
                          q= dim(Vt)[2],
                          lambda_prior = list(mean = rep(0, q),
                                              precision = rep(1, q))
                          ){
  n <- dim(X_batch)[1]; d <- dim(X_batch)[2]; q = dim(Vt)[2]
  Vt <- matrix(Vt, nrow= d, ncol =q)
  L_mle = batch_mle(X_batch, Vt, center, factor_family, q, weights)

  grad_v = mapply(PosSample_GradRowV,
                L_row = asplit(L_mle, 1),
                X_row = asplit(X_batch, 1),
                weight_row = asplit(weights, 1),
                MoreArgs = list(Vt = Vt, factor_family = factor_family,
                                dispersion = dispersion, sample_size = sample_size,
                                center = center, lambda_prior = lambda_prior))
  grad_v = matrix(rowSums(grad_v), nrow = d)

  grad_center = mapply(PosSample_GradRowCenter,
                L_row = asplit(L_mle, 1),
                X_row = asplit(X_batch, 1),
                weight_row = asplit(weights, 1),
                MoreArgs = list(Vt = Vt, factor_family = factor_family,
                                dispersion = dispersion, sample_size = sample_size,
                                center = center, lambda_prior = lambda_prior))
  grad_center = matrix(rowSums(grad_center), ncol = d)


  grad_dispersion = mapply(PosSample_GradRowPhi,
                       L_row = asplit(L_mle, 1),
                       X_row = asplit(X_batch, 1),
                       weight_row = asplit(weights, 1),
                       MoreArgs = list(Vt = Vt, factor_family = factor_family,
                                       dispersion = dispersion, sample_size = sample_size,
                                       center = center, lambda_prior = lambda_prior))

  grad_dispersion = matrix(rowSums(grad_dispersion), ncol = d)

  return (list(grad_V = grad_v, grad_center = grad_center,
               grad_dispersion = grad_dispersion))
}


# adam_update <- function(v, adam_control)
# v_dcenter = adam_control$beta1 * v_dcenter + (1 - adam_control$beta1) * grad$grad_center
# s_dcenter = adam_control$beta2 * s_dcenter + (1 - adam_control$beta2) * grad$grad_center^2
# vhat_dcenter = v_dcenter / (1 - adam_control$beta1 ^ adam_t)
# shat_dcenter = s_dcenter / (1 - adam_control$beta2 ^ adam_t)


# [function for cov experiments]
generate_cov <- function(n, d, L_prior, V_prior,center,
                         family, weights, phi){
  V0 = rmvnorm(d, mean = V_prior$mean,
               sigma = diag(V_prior$sigma),
               checkSymmetry = TRUE)
  L0 <- rmvnorm(n, mean = L_prior$mean,
                sigma = diag(1/L_prior$precision),
                checkSymmetry = TRUE)
  LV_identified <- efm_identifyLV(L0, V0)
  eta_LV = tcrossprod(LV_identified$L, LV_identified$V)
  eta0 = sweep(eta_LV, 2, center, '+')
  X = generate_data(family= family, eta = eta_LV,
                    dispersion = phi, weights =weights)
  list(X = X, L0 = L0, V0 = V0, phi = phi, center = center, weights = weights )
}

efm_identifyLV <- function(L, V){
  s <- sign(V[1, ])
  V_ <- sweep(V, 2, s, `*`)
  L_ <- sweep(L, 2, s, `*`)
  list(L = L_, V = V_)
}

':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL))
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL))
}


adam_update <- function(adam_param, grad_, adam_control,
                        lr_schedule, adam_t){

  v_candidate<- adam_control$beta1 * adam_param$v  + (1 - adam_control$beta1) * grad_
  s_candidate <- adam_control$beta2 * adam_param$s  + (1 - adam_control$beta2) * grad_^2

  #browser()
  if (!is.null(adam_param$ams_grad)){
    adam_param$ams_grad <- pmax(adam_param$ams_grad, s_candidate)
    s_update <- adam_param$ams_grad
  }else{
    s_update <-adam_param$s
  }

  adam_param$v <- v_candidate
  adam_param$s <- s_candidate

  vhat_dv = adam_param$v / (1 - adam_control$beta1^adam_t)
  shat_dv = s_update / (1 - adam_control$beta2^adam_t)

  param_update <-  lr_schedule * vhat_dv / (sqrt(shat_dv) + adam_control$epsilon)
  return(list(adam_param = adam_param, param_update = param_update))
}

