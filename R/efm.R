# [ Exponential Factor Model ]
#' @importFrom MASS ginv rnegbin
#' @importFrom glmnet glmnet
#' @importFrom stats coef dbinom dgamma dnbinom dnorm dpois gaussian rgamma rnorm rpois var rbinom
#' @importFrom matrixStats logSumExp rowLogSumExps
NULL
#> NULL
#'
# TODO:
# `add more annotation to utils.R
# `delete dependence on package mvtnorm, MASS, matrixStats
# `think about better stopping criteria


#' Generate data according to exponential family distribution given parameter `\eta`
#'
#' @param family Family object to specify factor loss.
#' @param eta `eta` to simulate data
#' @param weights Entry-wise weight.
#' @param dispersion dispersion parameter, effective for dispersed family, should be 1 for non-dispersed family.
#' @return simulated data from expoential family with mean `\mu = linkinv(\eta)`
#' @export
generate_data <-
  function (family,
            eta,
            weights = 1,
            dispersion = 1) {
    y_len <- prod(dim(eta))
    mu <- family$linkinv(eta)
    if (grepl("Negative Binomial", family$family)) family$family <- "negative.binomial"
    y <- switch(
      family$family,
      gaussian = rnorm(y_len, mu, sqrt(dispersion)),
      poisson = rpois(y_len, mu),
      binomial = rbinom(y_len, weights, mu),
      Gamma = rgamma(y_len , shape =  dispersion, scale = mu /  dispersion),
      negative.binomial = rnegbin(y_len, mu = mu, theta =  dispersion),
      stop("family `", family$family, "` not recognized")
    )
    dim(y) <- dim(eta)
    y
  }


family_initialize <- function (x, weights, family = gaussian()) {
  nr <- nrow(x); nc <- ncol(x)
  mustart <- NULL
  y <- c(x); nobs <- length(y) #; weights <- rep(1, nobs)
  eval(family$initialize)
  matrix(mustart, nrow = nr, ncol = nc)
}

#' Identify projection matrix `V` accrodng to `UD`
#'
#' @param V, Projection matrix to be identified
#' @export
efm_identify <- function(V){
    svd_efm = svd(V)
    svd_efm$u = efm_sign(svd_efm$u)
    return(sweep(svd_efm$u, 2, svd_efm$d, '*'))
}



batch_mle <- function(X_batch, Vt, center, factor_family, q, weights = 1) {
  n_sub = dim(X_batch)[1]
  d = dim(X_batch)[2]
  mu <- family_initialize(X_batch, weights = weights, factor_family)
  eta <- factor_family$linkfun(mu)
  L_mle <- eta[, 1:q, drop = FALSE]
  mu_eta <- matrix(factor_family$mu.eta(eta), ncol = d)
  var <- matrix(factor_family$variance(mu), ncol = d)
  is_inf_mu <- is.infinite(mu)
  S <- mu_eta  / var * mu_eta * weights
  S[is_inf_mu] <- 0
  Z <- sweep(eta, 2, center, '-' ) + (X_batch - mu) / mu_eta # working residuals
  Z[is_inf_mu] <- eta[is_inf_mu]
  for (i in 1:n_sub) {
    si <- sqrt(S[i, ])
    Vi <- sweep(Vt, 1, si, `*`)[, , drop = FALSE]
    L_mle[i, ] <-
      gsym_solve(crossprod(Vi), crossprod(Vi, (Z[i, ] * si)))
  }
  return(L_mle)
}


#' Using Monte Carlo to evaluate the marginalized exponential factor likelihood.
#'
#' @param Vt, Projection matrix for evaluation.
#' @param factor_family Family object to specify factor loss.
#' @param X Data to evaluate the likelihood.
#' @param weights Entry-wise weight.
#' @param sample_size number of Monte Carlo samples to estimate the gradients.
#' @param dispersion dispersion parameter, effective for dispersed family, should be 1 for non-dispersed family.
#' @param print_likeli boolen, to print out the simulated likelihood result
#' @param log_ boolen, to log the probability
#' @export
SML_neglikeli <- function(Vt,
                          factor_family,
                          X,
                          sample_size,
                          center = 0,
                          dispersion = 1,
                          weights = 1,
                          print_likeli = TRUE,
                          log_ = TRUE) {
  n <- dim(X)[1]; p <- dim(X)[2]; q <- dim(Vt)[2]
  Vt <- matrix(Vt, nrow = p, ncol = q)
  likeli_simu  <- matrix(0, nrow = n, ncol = sample_size)
  family_pdf <-
    pdf_calc(
      factor_family,
      weights = weights,
      dispersion = dispersion,
      log_ = log_
    )

  for (b in 1:sample_size) {
    L_sample <-  matrix(rnorm(n * q), nrow = n)
    eta_simu <- tcrossprod(L_sample, Vt)
    eta_simu = sweep(eta_simu, 2, center, "+")
    mu_simu <- factor_family$linkinv(eta_simu)
    likeli_simu [, b] <- rowSums(family_pdf(X,
                                            mu = mu_simu,
                                            weights = weights), na.rm =  TRUE)
  }
  neg_likeli = -sum(rowLogSumExps(likeli_simu) - log(sample_size))
  if (print_likeli) {
    print(neg_likeli)
  }
  return(neg_likeli)
}


#' Perform sign identification of SVD.
#' @param V  projection matrix of factor model
#' @export
efm_sign <- function (V) {
  s <- sign(V[1, ])
  return(V <- sweep(V, 2, s, `*`))
}

#' Adam optimization parameters
# the most important parameters are step_size, max_epoch and batch_size
#'
#' @param max_epoch number of optimization epoch.
#' @param batch_size number of data(rows) in each mini batch
#' @param step_size step size or learning rate for SGD.
#' @param beta1 controls the exponential decay rate used to scale the biased first moment estimate.
#' @param beta2 controls the exponential decay rate used to scale the biased second raw moment estimate.
#' @param rho learning rate decay through `step_size/(1 + 0.1 * t^{\rho})` with t being the iteration.
#' @param epsilon smoothing term to avoid division by zero.
#' @param abs_tol positive convergence tolerance `\epsilon`; the iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < `\epsilon`.
#' @return list of adam optimization parameters
#' @export
adam.control <- function(max_epoch = 10,
                         batch_size = 64,
                         step_size = 0.1,
                         rho = 1,
                         abs_tol = 1e-6,
                         beta1 = 0.9,
                         beta2 = 0.999,
                         epsilon = 10 ^ -8) {
  list(
    max_epoch = max_epoch,
    batch_size = batch_size,
    step_size = step_size,
    rho = rho,
    abs_tol = abs_tol,
    beta1 = beta1,
    beta2 = beta2,
    epsilon = epsilon
  )
}

#' EFM sampling parameters
#'
#' @param sample_size number of Monte Carlo samples to estimate the gradients
#' @param eval_size number of Monte Carlo samples to evaluate the marginal likelihood
#' @return list of efm sampling parameters
#' @export
sample.control <- function(sample_size = 50,
                           eval_size = 50) {
  list(sample_size = sample_size, eval_size = eval_size)
}


# NOTE: weight[i, j] == 0 means "don't use (i,j)"; is.na(x[i, j]) means "missing, estimate it"
#' Perform exponential factor model estimation.
#'
#' @param x Data to conduct factor inference
#' @param factor_family Family object to specify factor loss (see \code{family})
#' @param rank Rank of factor model.
#' @param weights Entry-wise weight.
#' @param algo Optimization algorithm to be chosen from 1). sml (simulated maximum likelihood) 2). ps (posterior sampling) 3). lapl (laplacian approxmation).
#' @param start matrix of initial projection matrix `V`
#' @param dispersion dispersion parameter, effective for dispersed family, should be 1 for non-dispersed family.
#' @param adam_control Adam optimization control parameters, (see \code{adam.control}).
#' @param sample_control EFM sampling control parameters, (see \code{sample.control})
#' @param eval_likeli boolen, set to false to skip marginal likelihood evaluation
#' @param identify_ boolen, set to true to enable identifiability at every iteration
#' @return Estimation of projection matrix `V`.
#' @export
efm <- function(x,
                factor_family,
                rank,
                weights = 1,
                algo = 'lapl',
                start = NULL,
                dispersion = 1,
                center = matrix(0, ncol = ncol(x)),
                adam_control = adam.control(
                  max_epoch = 5,
                  batch_size = 32,
                  step_size = 0.1,
                  rho = 0,
                  abs_tol = 1e-6,
                  beta1 = 0.9,
                  beta2 = 0.999,
                  epsilon = 10 ^ -8),# typo
                sample_control = sample.control(sample_size = 50, eval_size = 50),
                eval_likeli = FALSE,
                identify_ = FALSE) {
  n <- nrow(x); p <- ncol(x)
  rank <- as.integer(rank)
  if (rank <= 0)
    stop("invalid non-positive rank")
  if (length(weights) == 1)
    weights <- rep(weights, n * p)
  else {
    if (is.vector(weights)) {
      if (length(weights) != n)
        stop("inconsistent number of weights")
      else
        weights <- rep(weights, p)
    } else {
      if (nrow(weights) != n && ncol(weights) != p)
        stop("inconsistent number of weights")
    }
  }
  dim(weights) <- dim(x)

  if (is.null(start)) {
    # [ initialize through SVD]
    mu <- family_initialize(x, weights, factor_family)
    eta <- factor_family$linkfun(mu)
    se <- svd(eta, nu = rank, nv = rank)
    Vt <- sweep(se$v, 2, se$d[1:rank], `*`)

  } else {
    Vt = start
    if (nrow(Vt) != p || ncol(Vt) != rank)
      stop("dimensions of V are inconsistent")
  }

  v_dv = matrix(0, nrow = rank, ncol = p); s_dv = matrix(0, nrow = rank, ncol = p)
  v_dcenter = matrix(0, ncol = p); s_dcenter = matrix(0, ncol = p)
  like_list <- rep(0, adam_control$max_epoch * as.integer(n / adam_control$batch_size))

  for (epoch in 1:adam_control$max_epoch) {
    for (k in 1:as.integer(n / adam_control$batch_size)) {
      sample_index = sample(1:n, adam_control$batch_size, replace = TRUE)
      x_batch <- x[sample_index,]
      weight_batch <- weights[sample_index,]

      # [Finding gradient moment using posterior Lambda]
      grad <- switch(
        algo,
        lapl = Lapl_grad(X_batch = x_batch, Vt = Vt, factor_family = factor_family, weights = weight_batch, dispersion = dispersion, center = center),
        sml = SML_grad(
          Vt = Vt,
          factor_family = factor_family,
          X = x_batch,
          q = rank,
          center = center,
          sample_size = sample_control$sample_size,
          dispersion = dispersion,
          weights = weight_batch
        ),
        fastgaussian = FastGaussian_grad(x_batch, Vt, factor_family, dispersion, weight_batch),
        ps = PosSample_grad(
          X_batch = x_batch,
          Vt = Vt,
          center = center,
          factor_family = factor_family,
          dispersion = dispersion,
          weights = weight_batch,
          sample_size = sample_control$sample_size
        ),
        stop("Algo `", algo, "` not recognized")
      )

      grad = lapply(grad, "*" ,n / adam_control$batch_size)

      # [Adam lr decay]
      adam_t = (epoch - 1) * as.integer(n / adam_control$batch_size) + k
      lr_schedule = adam_control$step_size / (1 + 0.1 * adam_t ^ adam_control$rho)

      # [Adam V step]
      v_dv = adam_control$beta1 * v_dv + (1 - adam_control$beta1) * grad$grad_V
      s_dv = adam_control$beta2 * s_dv + (1 - adam_control$beta2) * grad$grad_V^2
      vhat_dv = v_dv / (1 - adam_control$beta1 ^ adam_t)
      shat_dv = s_dv / (1 - adam_control$beta2 ^ adam_t)

      # [ Adam center step ]
      v_dcenter = adam_control$beta1 * v_dcenter + (1 - adam_control$beta1) * grad$grad_center
      s_dcenter = adam_control$beta2 * s_dcenter + (1 - adam_control$beta2) * grad$grad_center^2
      vhat_dcenter = v_dcenter / (1 - adam_control$beta1 ^ adam_t)
      shat_dcenter = s_dcenter / (1 - adam_control$beta2 ^ adam_t)

      Vt_update =  lr_schedule * t(vhat_dv / (sqrt(shat_dv) + adam_control$epsilon))
      center_update = lr_schedule * vhat_dcenter / (sqrt(shat_dcenter) + adam_control$epsilon)

      Vt = Vt - Vt_update
      center = center - center_update
      # print(center_update)
      #print(center)

      #[Identifiability]
      if (identify_) {
        svd_efm = svd(Vt)
        svd_efm$u = efm_sign(svd_efm$u)
        Vt = sweep(svd_efm$u, 2, svd_efm$d, '*')
      }

      # [evaluate marginal likelihood]
      if (eval_likeli) {
        like_list[adam_t] = SML_neglikeli(
          Vt,
          factor_family,
          x,
          center = center,
          sample_control$eval_size,
          dispersion = dispersion,
          weights = weights,
          print_likeli = TRUE
        )
        plot(like_list[1:adam_t])
      }

      # [early stopping]
      update_norm = mean((Vt_update) ^ 2)
      #print(update_norm)
      stop_flag = update_norm < adam_control$abs_tol
      if ((epoch >= 2) && (stop_flag)) {
        print('find optimal points')
        return(list(
          V = Vt,
          family = factor_family,
          efm_it = adam_t,
          center = center,
          like_list = like_list[1:adam_t]
        ))
      }
    } # end of one epoch
  } # end of max epoch
  return(list(
    V = Vt,
    family = factor_family,
    center = center,
    efm_it = adam_t,
    like_list = like_list[1:adam_t]
  ))

} # end of function
