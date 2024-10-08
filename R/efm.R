# [ Exponential Factor Model ]
#' @importFrom MASS ginv rnegbin
#' @importFrom stats coef dbinom dgamma dnbinom dnorm dpois gaussian rgamma rnorm rpois var rbinom
#' @importFrom matrixStats logSumExp rowLogSumExps
NULL
#> NULL
#'
# TODO:
# `add more annotation to utils.R
# `delete dependence on package mvtnorm, MASS, matrixStats
# `think about better stopping criteria


# TODO: how to generate when dispersion<1?
rqpoisson <- function(n, mu, dispersion){rnbinom(n = n, mu = mu, size = mu/(dispersion-1))}
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
    n = dim(eta)[1]; p = dim(eta)[2]
    y_len <- prod(dim(eta))
    mu <- family$linkinv(eta)
    if (grepl("Negative Binomial", family$family)) family$family <- "negative.binomial"
    #rnorm(100, mean = 1:100, sd = c(rep(0.001,99), 100)),
      # dispersion with same dim on num_generate, mu, phi, works
    # dispersion <- switch (length(dispersion),
    #   '1' = matrix(dispersion, nrow = n, ncol = p),
    #   'p' = do.call("rbind", replicate(n, dispersion, simplify = FALSE)),
    #   'n' = do.call("cbind", replicate(p, dispersion, simplify = FALSE)),
    #   y_len = dispersion
    # )
    dispersion_len <- length(dispersion)
    if (dispersion_len == 1){
      dispersion <- matrix(dispersion, nrow = n, ncol = p)
    }else if (dispersion_len == p){
      dispersion <-do.call("rbind", replicate(n, dispersion, simplify = FALSE))
    }else if (dispersion_len == n){
      dispersion <-do.call("cbind", replicate(p, dispersion, simplify = FALSE))
    }

    y <- switch(
      family$family,
      gaussian = rnorm(y_len, mu, sqrt(dispersion)),
      poisson = rpois(y_len, mu),
      quasipoisson = rqpoisson(y_len, mu, dispersion),
      binomial = rbinom(y_len, weights, mu),
      Gamma = rgamma(y_len , shape =  1/dispersion, scale = mu * dispersion), # V(x) = dispersion * mean(x)^2
      negative.binomial = rnegbin(y_len, mu = mu, theta = dispersion),
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
      gsym_solve(symm_eigen(crossprod(Vi)), crossprod(Vi, (Z[i, ] * si)))
  }
  return(L_mle)
}

#' @export
#' Var(X) = E_L[Var(X|L) + Var(E(X|L))]
marg_var <- function (center, V, family, ngq = 15, dispersion =1,
                      L_prior = list(mean = rep(0, q), precision = rep(1, q)), ...) {
  if (!is.vector(center)) center = c(center)
  gq <- gaussquadr::GaussQuad$new("gaussian", ngq, ncol(V), ...)
  gi <- gq$clone()$location_scale(L_prior$mean, 1/L_prior$precision, squared = TRUE)


  fam <- if (is.function(family)) family() else family
  ev <- gi$E(function (l) fam$variance(fam$linkinv(center + V %*% l)) * dispersion )
  mv <- gi$E(function (l) fam$linkinv(center + V %*% l))
  ve <- gi$E(function (l) tcrossprod(fam$linkinv(center + V %*% l) - mv))
  matrix(ve, nrow = nrow(V)) + diag(ev)
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
                          L_prior = list(mean = rep(0, q),
                                         precision = rep(1, q)),
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
      family = factor_family,
      dispersion = dispersion,
      weights = weights,
      log_ = log_
    )
  dispersion_len <- length(dispersion)
  if (dispersion_len == 1){
    dispersion_eval <- matrix(dispersion, nrow = n, ncol = p)
  }else if (dispersion_len == p){
    dispersion_eval <-do.call("rbind", replicate(n, dispersion, simplify = FALSE))
  }else if (dispersion_len == n){
    dispersion_eval <-do.call("cbind", replicate(p, dispersion, simplify = FALSE))
  }

  for (b in 1:sample_size) {
    L_sample <-  rmvnorm(n, mean = L_prior$mean,
                         sigma = diag(1/L_prior$precision),
                         checkSymmetry = TRUE)

    eta_simu <- tcrossprod(L_sample, Vt)
    eta_simu = sweep(eta_simu, 2, center, "+")
    mu_simu <- factor_family$linkinv(eta_simu)

    family_density <- family_pdf(X,
               mu = mu_simu,
               weights = weights,
               dispersion = dispersion_eval)
    family_density[is.na(family_density)] <- max(family_density[!is.na(family_density)])

    likeli_simu [, b] <- rowSums(family_density)
  }
  neg_likeli = -sum(rowLogSumExps(likeli_simu) - log(sample_size))
  #neg_likeli = -sum(likeli_simu)/sample_size
  if (print_likeli) {
    print(paste('loss(negative quasi-likelihood) is', neg_likeli))
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


#' @param X Data to conduct factor inference
#' @param q Rank of factor model.
#' @param family Family object to specify factor loss (see \code{family})
#' @param ngq Number of nodes used in gaussian quadrature
#' @param weights Entry-wise weight.
#' @param lambda_prior Prior on Lambda, a list containing mean and precision
#' @param control Optimization configuration, same to glm.control, (see \code{glm.control})
#' @param eval_likeli Boolen, set to false to skip marginal likelihood evaluation
#' @param eval_size MC sample size used to evaluate negative likelihood
#' @param start Initialization point, a list containing V, Phi and center
#' @return Estimation of projection matrix `V`.
fa_gqem <- function (X, q, ngq, family = gaussian(), weights,
                      lambda_prior = list(mean = rep(0, q), precision = rep(1, q)),
                      control = list(), eval_size = 500,
                      eval_likeli = FALSE, start = NULL, ...) {
  n <- nrow(X); p <- ncol(X)
  control <- do.call("glm.control", control)
  fam <- if (is.function(family)) family() else family

  # [initialize weights]
  if (missing(weights)) weights <- 1
  if (length(weights) == 1)
    weights <- rep(weights, n * p)
  else {
    if (is.vector(weights)) {
      if (length(weights) != p)
        stop("inconsistent number of weights")
      else
        weights <- weights[gl(p, n)] # kronecker(weights, rep(1, n))
    } else {
      if (nrow(weights) != n && ncol(weights) != p)
        stop("inconsistent number of weights")
    }
  }
  dim(weights) <- dim(X)

  # [initialize iteration]
  phi_flag <- check_DispersionUpdate(fam)
  if (is.null(start)) {
    mu <- family_initialize(X, weights, fam)
    eta <- fam$linkfun(mu)
    scale_eta <- scale(eta, scale = FALSE) # center
    alpha <- attr(scale_eta, "scaled:center")
    S_ <- cov(scale_eta)
    ss_ <- symm_eigen(S_)
    V <- sweep(ss_$vectors[, 1:q, drop = FALSE], 2, sqrt(ss_$values[1:q]), `*`)
    if (phi_flag) Phi <- apply(weights * (X - mu) ^ 2 / fam$variance(mu), 2, mean)
  } else {
    V <- start$Vt
    alpha <- c(start$center)
    if (phi_flag) {
      Phi <- c(start$dispersion)
      if (length(Phi) == 1) Phi <- rep(Phi, p)
    }
  }
  if (!phi_flag) Phi <- rep(1, p)
  Phi[Phi <= 0] <- 0.01

  # [initialize integration]
  gq <- gaussquadr::GaussQuad$new("gaussian", ngq, q, ...)
  m <- length(gq$weights)
  lw_ref <- -.5 * apply(gq$nodes, 1, norm2) ^ 2 # approx posterior
  dev <- Inf

  # [initialize mc likelihood]
  like_list <- rep(0, control$maxit + 1)

  if (eval_likeli) {
    like_list[1] <- SML_neglikeli(V = V, factor_family = fam, X = X,
                                sample_size = eval_size, L_prior = lambda_prior,
                                center = alpha, dispersion = Phi, weights = weights)}
  eval_time <- 0
  for (it in 1:control$maxit) {
    # [ E-step: lambda | x ]
    H <- matrix(0, nrow = p, ncol = (q + 1) ^ 2)
    z <- matrix(0, nrow = p, ncol = q + 1)
    res <- numeric(p)
    dev_new <- 0
    for (i in seq_len(n)) {
      bs <- bsglm(V, X[i, ], lambda_prior, offset = alpha, family = fam,
                  weights = weights[i, ] / Phi, return_hessian = TRUE,
                  control = list(maxit = 1))

      bs$hessian$values <- 1 / bs$hessian$values # invert first
      if (!all(bs$hessian$valid))
        warning("only ", sum(bs$hessian$valid), " valid eigenvalues")
      # bs$hessian$values[bs$hessian$values <=0] <- min(bs$hessian$values[bs$hessian$values >0])
      # bs$hessian$values[is.infinite(bs$hessian$values)] <- max(bs$hessian$values[!is.infinite(bs$hessian$values)])

      gi <- gq$clone()$location_scale(bs$coef, bs$hessian, squared = TRUE)
      # refine weights with importance
      mu <- fam$linkinv(tcrossprod(V, gi$nodes) + alpha)
      devi <- apply(mu, 2, function (mui)
                    sum(fam$dev.resids(X[i,], mui, weights[i, ] / Phi)))
      devi <- devi + apply(gi$nodes, 1, norm2) ^ 2 # add log prior
      gi$reweight(-.5 * devi - lw_ref, log_weights = TRUE, normalize = TRUE)

      for (j in seq_len(p)) { # fit column regressions
        etaj <- drop(gi$nodes %*% V[j, ] + alpha[j])
        muj <- fam$linkinv(etaj)
        varmuj <- fam$variance(muj)
        wj <- weights[i, j] * gi$weights
        Wj <- wj * varmuj
        rj <- wj * (X[i, j] - muj)
        if (!check_canonical(fam)) { # adjust weights?
          thetaetaj <- fam$mu.eta(etaj) / varmuj
          Wj <- Wj * thetaetaj ^ 2
          rj <- rj * thetaetaj
        }
        xj <- cbind(1, gi$nodes)
        z[j, ] <- z[j, ] + crossprod(xj, Wj * etaj + rj)
        #if(sum(is.na(z))!=0){browser()}
        H[j, ] <- H[j, ] + c(crossprod(xj, xj * Wj))
        res[j] <- res[j] + sum(fam$dev.resids(rep(X[i, j], m), muj, wj))
      }
    }

    # [ M-step: fit alpha, V, and Phi ]
    for (j in seq_len(p)) {
      betaj <- gsym_solve(symm_eigen(matrix(H[j, ], nrow = q + 1)), z[j, ])
      alpha[j] <- betaj[1]; V[j, ] <- betaj[-1]
    }
    if (phi_flag) Phi <- res / n
    dev_new <- sum(res / Phi + n * log(Phi))

    # [ update mc likelihood]
    if (eval_likeli) {
      start <- Sys.time()
      like_list[it + 1] <- SML_neglikeli(V = V, factor_family = fam, X = X,
                                sample_size = eval_size, L_prior = lambda_prior,
                                center = alpha, dispersion = Phi, weights = weights)
      plot(like_list[1:(it +1)])
      end <- Sys.time()
      eval_time <- eval_time + end -start
    }

    if (control$trace) {
      message("[", it, "] dev = ", dev_new)
    }
    # if (it > 1 && abs((dev_new - dev) / (dev + .1)) < control$epsilon) break
    # dev <- dev_new
  }
  sv <- svd(V, nv = 0); V <- sweep(sv$u, 2, sv$d * sign(sv$u[1,]), "*")
  list(center = alpha, V = V, like_list = like_list[like_list!=0], eval_time = eval_time,
       family = fam, ngq = ngq, deviance = dev, dispersion = Phi)
}

# NOTE: weight[i, j] == 0 means "don't use (i,j)"; is.na(x[i, j]) means "missing, estimate it"
#' Perform exponential factor model estimation.
#'
#' @param x Data to conduct factor inference
#' @param factor_family Family object to specify factor loss (see \code{family})
#' @param rank Rank of factor model.
#' @param weights Entry-wise weight.
#' @param algo Optimization algorithm to be chosen from 1). sml (simulated maximum likelihood) 2). ps (posterior sampling) 3). lapl (laplacian approxmation).
#' @param start Initialization point, a list containing V, Phi and center
#' @param lambda_prior Prior on Lambda, a list containing mean and precision
#' @param adam_control Adam optimization control parameters, (see \code{adam.control}).
#' @param sample_control EFM sampling control parameters, (see \code{sample.control})
#' @param em_control Same to glm.control, used for fa_gqem, (see \code{glm.control})
#' @param ngq Number of nodes used in gaussian quadrature
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
                lambda_prior = list(mean = rep(0, q),
                                    precision = rep(1, q)),
                adam_control = adam.control(
                  max_epoch = 5,
                  batch_size = 32,
                  step_size = 0.1,
                  rho = 0,
                  abs_tol = 1e-6,
                  beta1 = 0.9,
                  beta2 = 0.999,
                  epsilon = 10^-8),
                sample_control = sample.control(sample_size = 50, eval_size = 500),
                em_control = list(),
                ngq = 15,
                eval_likeli = FALSE,
                identify_ = FALSE) {

  rank <- as.integer(rank)
  if (rank <= 0)
    stop("invalid non-positive rank")
  if (algo == 'em'){
    result <- fa_gqem(X = x, q = rank, ngq, family = factor_family, weights = weights,
                      lambda_prior = lambda_prior,
                      control = em_control,  eval_size = sample_control$eval_size,
                      eval_likeli = eval_likeli, start = start)
    return(result)
  }


  n <- nrow(x); p <- ncol(x)
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
    scale_eta <- scale(eta, scale = FALSE) # center
    center <- attr(scale_eta, "scaled:center")
    S_ <- cov(eta)
    ss_ <- symm_eigen(S_)
    Vt <- sweep(ss_$vectors[, 1:q, drop = FALSE], 2, sqrt(ss_$values[1:q]), `*`)
    dispersion = apply(weights * (x - mu)^2/factor_family$variance(mu), 2, mean)
  } else {
    Vt = start$Vt
    center = start$center
    dispersion = start$dispersion
    if (length(dispersion) ==1)(dispersion = rep(dispersion, p))
    if (nrow(Vt) != p || ncol(Vt) != rank)
      stop("dimensions of V are inconsistent")
  }
  phi_flag = check_DispersionUpdate(factor_family)
  dispersion[dispersion<=0] = 0.1
  if(phi_flag == FALSE){dispersion = 1}



  total_iter <- adam_control$max_epoch * as.integer(n / adam_control$batch_size)
  like_list <- rep(0, total_iter + 1)

  like_list[1] <-SML_neglikeli(V = Vt, factor_family = factor_family, X = x, sample_size = sample_control$eval_size, L_prior = lambda_prior,
                                center = center, dispersion = dispersion, weights = weights)

  eval_time <-0

  adam_V <- list(v = matrix(0, nrow = p, ncol = rank),
                s = matrix(0, nrow = p, ncol = rank),
                ams_grad = matrix(0, nrow = p, ncol = rank))
  adam_center <- list(v = matrix(0, ncol = p),
                     s = matrix(0, ncol = p),
                     ams_grad = matrix(0, ncol = p))
  adam_phi <- list(v = matrix(0, ncol = p),
                   s = matrix(0, ncol = p),
                   ams_grad =  matrix(0, ncol = p))

  for (epoch in 1:adam_control$max_epoch) {
    for (k in 1:as.integer(n / adam_control$batch_size)) {
      sample_index = sample(1:n, adam_control$batch_size, replace = TRUE)
      x_batch <- x[sample_index,]
      weight_batch <- weights[sample_index,]

      # [Finding gradient moment using posterior Lambda]
      grad <- switch(
        algo,
        lapl = Lapl_grad(X_batch = x_batch,
                         Vt = Vt, factor_family = factor_family,
                         weights = weight_batch,
                         dispersion = dispersion,
                         center = center,
                         lambda_prior = lambda_prior),
        sml = SML_grad(
          Vt = Vt, factor_family = factor_family,
          X = x_batch, q = rank,center = center,
          sample_size = sample_control$sample_size,
          dispersion = dispersion, weights = weight_batch,
          lambda_prior = lambda_prior
        ),
        ps = PosSample_grad(
          X_batch = x_batch, Vt = Vt,
          center = center, factor_family = factor_family,
          dispersion = dispersion, weights = weight_batch,
          sample_size = sample_control$sample_size,
          lambda_prior = lambda_prior
        ),
        stop("Algo `", algo, "` not recognized")
      )

      grad = lapply(grad, "*" ,n / adam_control$batch_size)
      # [Adam lr decay]
      adam_t = (epoch - 1) * as.integer(n / adam_control$batch_size) + k
      lr_schedule = adam_control$step_size / (1 + 0.1 * adam_t ^ adam_control$rho)

      # [Adam V step]
      c(adam_V, Vt_update) := adam_update(adam_V, grad$grad_V, adam_control, lr_schedule, adam_t)
      c(adam_center, center_update) := adam_update(adam_center, grad$grad_center, adam_control, lr_schedule, adam_t)
      c(adam_phi, dispersion_update) := adam_update(adam_phi, grad$grad_dispersion, adam_control, lr_schedule, adam_t)

      # [Update steps]
      Vt = Vt - Vt_update
      center = center - center_update


      if (phi_flag == TRUE){
        dispersion = dispersion - c(dispersion_update)
        invalid_flag = dispersion<0
        dispersion[invalid_flag] = 0.1}


      #[Identifiability]
      if (identify_) {
        svd_efm = svd(Vt)
        svd_efm$u = efm_sign(svd_efm$u)
        Vt = sweep(svd_efm$u, 2, svd_efm$d, '*')}

      # [evaluate marginal likelihood]
      if (eval_likeli) {
        start <- Sys.time()
        like_list[adam_t + 1] <- SML_neglikeli(V = Vt, factor_family = factor_family, X = x, sample_size = sample_control$eval_size, L_prior = L_prior,
                                center = center, dispersion = dispersion, weights = weights)
        plot(like_list[1: (adam_t + 1)])
        end <- Sys.time()
        eval_time <- eval_time + end - start
      }

      # [early stopping, setting update_norm = Inf to disable]
      update_norm = Inf
      stop_flag = update_norm < adam_control$abs_tol
      if ((epoch >= 2) && (stop_flag)) {
        print('find optimal points')
        return(list(
          V = Vt, family = factor_family,
          efm_it = adam_t, center = center,
          dispersion = dispersion, like_list = like_list[1:adam_t],
          algo = algo))
      }
    } # end of one epoch
  } # end of max epoch
  return(list(
    V = Vt, family = factor_family,
    center = center, efm_it = adam_t,
    dispersion = dispersion, like_list = like_list,
    algo = algo, eval_time = eval_time
  ))

} # end of function


