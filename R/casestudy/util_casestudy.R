
# max_epoch = 1;sample_size = 100; beta1 = 0.9; beta2 = 0.999;
# epislon = 10^-8; phi_star = phi_star; rho = 1; scale_weights = 1
# batch_opti<-function(Vt, batch_size, step_size, X, factor_family, q, max_epoch = 1,
#                      sample_size = 100, beta1 = 0.9, beta2 = 0.999,epislon = 10^-8,
#                       phi_star = 1, rho = 1, scale_weights = 1){
batch_opti<-function(Vt, batch_size, step_size, X, factor_family, q, max_epoch = 1,
                     sample_size = 100, beta1 = 0.9, beta2 = 0.999,epislon = 10^-8,
                     phi_star = 1, rho = 1, scale_weights = 1, like_threshold =0 ,sample_random = TRUE){
  # [Adam initialization]
  n = dim(X)[1]
  d = dim(X)[2]
  family_pdf = pdf_calc(factor_family, glm_weights = 1, phi= phi_star, log = TRUE)
  v_dv = matrix(0, nrow = q, ncol = d)
  s_dv = matrix(0, nrow = q, ncol = d)
  like_list = rep(0, max_epoch * as.integer(n/batch_size))
  #like_list2 = rep(0, max_epoch * as.integer(n/batch_size))
  for (epoch in 1:max_epoch){
    for(k in 1:as.integer(n/batch_size)){
      if(sample_random ==TRUE){
        sample_index = sample(1:n, batch_size)
      }else{
        batch_size = n
        sample_index = 1:n
      }
      X_batch = X[sample_index,]
      n_sub = dim(X_batch)[1]
      grad = matrix(0, nrow = q, ncol = d)


      # [Finding the stationary point of L]
      L_mle = batch_mle(X_batch, Vt, factor_family, q)
      mu_mle = factor_family$linkinv(tcrossprod(L_mle,Vt))


      # [Computing Laplacian sampling parameters]
      if (factor_family$family=='gaussian'){
        hessian = scale_weights/phi_star^2 * crossprod(Vt) # for gaussian
        mu_pos = tcrossprod(ginv(ginv(hessian) + diag(1, nrow= q)), L_mle)
        Sigma_pos = ginv(diag(1, nrow= q) + hessian)
        Sigma_chol = chol(Sigma_pos)}
      else{
        Weight_matrix = cbind(Weight_row = tcrossprod(L_mle, Vt), L_row = L_mle) # for non-gaussian
        rmv_norm = function(Weight_row) {
          L_row = matrix(Weight_row[(d+1):(d+q)], nrow= 1)
          Weight_row = Weight_row[1:d]
          mu_hat = factor_family$linkinv(Weight_row)
          diag_term = diag(scale_weights * factor_family$variance(mu_hat), nrow = d)
          hessian = tcrossprod(crossprod(Vt, diag_term), t(Vt))
          Sigma_pos = ginv(diag(1, nrow= q) + hessian) + diag(1e-08, nrow = q)
          Sigma_chol = chol(Sigma_pos)
          mu_pos = tcrossprod(ginv(ginv(hessian) + diag(1, nrow= q)), L_row)
          return(mu_pos + tcrossprod(Sigma_chol, matrix(rnorm(1 * q), nrow = 1)))
        }
      }


      like_temp = matrix(0, nrow = n, ncol = sample_size)
      for (b in 1:sample_size){
        #[gaussian sample]
        if (factor_family$family=='gaussian'){
          L_sample = t(mu_pos + tcrossprod(Sigma_chol, matrix(rnorm(n_sub * q), nrow = n_sub)))
          grad = grad + (scale_weights /phi_star^2 * crossprod(L_sample, factor_family$linkinv(tcrossprod(L_sample, Vt)) - X_batch))
          L_sample2 =  matrix(rnorm(n*q), nrow = n)
          like_temp[,b] = rowSums(family_pdf(X, mu = factor_family$linkinv(tcrossprod(L_sample2, Vt)), glm_weights =1))

        }else{
          #[poisson sample]
          L_sample = t(matrix(apply(Weight_matrix, 1,rmv_norm), ncol = batch_size))
          grad = grad + (scale_weights* crossprod(L_sample, factor_family$linkinv(tcrossprod(L_sample, Vt)) - X_batch))
          L_sample2 =  matrix(rnorm(n*q), nrow = n)
          like_temp[,b] = rowSums(family_pdf(X, mu = factor_family$linkinv(tcrossprod(L_sample2, Vt)), glm_weights = 1))
        }
      }
      grad = grad / sample_size  * n/batch_size
      #[ Adam step]
      v_dv = beta1 * v_dv + (1-beta1)*grad
      s_dv = beta2 * s_dv + (1-beta2)*grad^2
      adam_t = (epoch-1)*as.integer(n/batch_size) + k
      vhat_dv = v_dv/(1-beta1^adam_t)
      shat_dv = s_dv/(1-beta2^adam_t)

      Vt = Vt - step_size/(1 + 0.1 * adam_t^rho) * t(vhat_dv/(sqrt(shat_dv) + epislon))



      #[Identifiability]
      svd_efm = svd(Vt)
      negative_flag = c(1:q)[svd_efm$u[1,]<0]
      svd_efm$u[,negative_flag] = matrix(-svd_efm$u[,negative_flag])
      Vt = tcrossprod(svd_efm$u, diag(svd_efm$d))
      #print(Vt)


      #[computing the marginal likelihood]
      if (factor_family$family=='gaussian'){
        #like_list[adam_t] = -sum(dmvnorm(X, sigma = tcrossprod(Vt) + diag(phi_star^2,d), log = TRUE))
        like_list[adam_t] = -sum(rowLogSumExps(like_temp)- log(sample_size))
      }else{
        like_list[adam_t] = -sum(rowLogSumExps(like_temp)- log(sample_size))
      }
      #plot(like_list[1:adam_t])

      #points(like_list2[1:adam_t], col = 'red')
      #[earlier stopping]
      if (epoch>=2){
        stop_flag = abs(like_list[adam_t] - like_list[adam_t-1])/like_list[adam_t-1]<= like_threshold

        if (stop_flag ==TRUE){
          L_mle = batch_mle(X, Vt, factor_family, q)
          print('find optimal points')
          return(list(like_list = like_list, V = Vt, L = L_mle, family = factor_family))}
      }

    } #end of batch
  } # end of epoch
  L_mle = batch_mle(X, Vt, factor_family, q)
  L_pos = apply(L_mle,1, comput_mupos, Vt, factor_family = factor_family)
  return(list(like_list = like_list, V = Vt, L = t(L_pos), family = factor_family))
}
