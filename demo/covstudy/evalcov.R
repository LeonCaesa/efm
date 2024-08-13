library(efm)
#setwd('/projectnb2/dmfgrp/efm')
#source('../R/utils.R')
library(MASS)
library(tidyverse)
n_gq <- 15



names_algo <- c('lapl', 'sml', 'ps')
#
# metric_cov <- function(esti_cov, actual_cov){
#   d <- dim(esti_cov)[1]
#   product = esti_cov %*% ginv(actual_cov)
#   l2_frobenius <- sum((esti_cov - actual_cov)^2)
#   l1_entropy <- sum(diag(product)) - log (det(product)) - d
#   l2_normalized <- sum((diag(product) -1)^2)* d^(-1/2)
#   return(list(l2_frobenius, l1_entropy, l2_normalized))
# }

metric_cov <- function(esti_cov, actual_cov){
  d <- dim(esti_cov)[1]
  eigen_frobenius <- eigen(esti_cov - actual_cov, only.values = TRUE)
  l2_frobenius <- sqrt(sum(eigen_frobenius$values^2))

  product = esti_cov %*% ginv(actual_cov)
  l1_entropy <- sum(diag(product)) - log (det(product)) - d

  sigma_half <- chol(ginv(actual_cov))
  inside_norm <- sigma_half %*% (esti_cov - actual_cov) %*% t(sigma_half)
  eigen_normalized <- eigen(inside_norm, only.values = TRUE)
  l2_normalized = sqrt(sum( eigen_normalized$values^2) ) *d^(-1/2)
  return(list(l2_frobenius, l1_entropy, l2_normalized))
}


# TODO: check why l2 normalized >1
# load("/projectnb/dmfgrp/efm/CovResult1209/negbinom(20)/_total_466.RData")
# actual_cov <- true_cov
# #esti_cov <- fagqem_esti
# esti_cov <- naive_esti
#
# # [method 1]
# sigma_half <- chol(ginv(actual_cov))
# inside_norm <- sigma_half %*% (esti_cov - actual_cov) %*% t(sigma_half)
# eigen_normalized <- eigen(inside_norm, only.values = TRUE)
# l2_normalized = sqrt(sum( eigen_normalized$values^2) ) *d^(-1/2)
#
# # [method 1]
# input_diff <- esti_cov %*% ginv(actual_cov) - diag(1, d)
# sum(diag(input_diff)^2) * d^(-1/2)

#d_list = seq(10, 100, by = 10); n = 500
# factor_family = quasipoisson()
# load_dir = '/projectnb/dmfgrp/efm/CovResult1209/quasipoisson/'
# load_dir = '/projectnb/dmfgrp/efm/CovResult/quasipoisson/'


# factor_family = negative.binomial(2)
# load_dir = '/projectnb/dmfgrp/efm/CovResult/negbinom_phi2/'

# factor_family = negative.binomial(20)
# load_dir = '/projectnb/dmfgrp/efm/CovResult1209/Negative Binomial(20)/'
#load_dir = '/projectnb/dmfgrp/efm/CovResult1001/Negative Binomial(20)_Fagqem/Negative Binomial(20)/'

# factor_family = binomial()
# load_dir = '/projectnb/dmfgrp/efm/CovResult1209/binomial/'

# factor_family = poisson()
# load_dir = '/projectnb/dmfgrp/efm/CovResult1209/poisson/'
#load_dir = '/projectnb/dmfgrp/efm/CovResult1001/poisson/'

#d_list = seq(16, 1000, by = 100); n = 756
#d_list = seq(16, 716, by = 100); n = 756
d_list = seq(66, 466, by = 50); n = 756


# [To store Cov Error]
error_matrix <- matrix(ncol = 5)
neglikeli_matrix<-matrix(ncol = 5)


# [To store likelihood decrease]
adam_control <- adam.control(max_epoch = 20, batch_size = 128,
                             step_size = 0.2, rho = 0, abs_tol = 1e-8,
                             beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8)
total_steps <- adam_control$max_epoch * as.integer(n / adam_control$batch_size)
neglikeli_array <- array( dim = c(total_steps, ncol = 5, 5))
colnames(neglikeli_array) <- c('naive', 'fqgqem', names_algo)

metric_names = c('l2frobenius', 'l2entrophy', 'l2normalized')



family_list <- c('quasipoisson', 'negbinom(20)', 'poisson', 'binomial')

for (family_name in family_list){

    load_dir = paste('/projectnb/dmfgrp/efm/CovResult1209', family_name, '', sep = '/')


    factor_family <- switch(family_name, 'quasipoisson' = quasipoisson(),
                           'negbinom(20)' = negative.binomial(20),
                           'poisson' = poisson(),
                           'binomial' = binomial())

    # [To store Cov Error]
    error_matrix <- matrix(ncol = 5)
    neglikeli_matrix<-matrix(ncol = 5)

    for (repeat_idx in 1:20){

      d_count = 1
        for (d in d_list){

            print( c(family_name, repeat_idx, d))

            truth_name <- paste(paste(load_dir, 'truth', d, repeat_idx, sep = '_'), '.RData', sep = '')
            load(truth_name)
            fagqem_name <- paste(paste(load_dir, 'fageqm', d, repeat_idx, sep = '_'), '.RData', sep = '')
            load(fagqem_name)

            if (grepl("Negative Binomial", factor_family$family)) truth$phi <- 1

            true_cov <- marg_var(truth$center, truth$V, family = factor_family,
                                 ngq = 15, dispersion = truth$phi, L_prior = truth$L_prior)
            fagqem_esti <- marg_var(efm_fagqem$center, efm_fagqem$V,
                                    family = factor_family, ngq = 15,
                                    dispersion = efm_fagqem$dispersion, L_prior = truth$L_prior)
            naive_esti <- cov(truth$X/truth$weights)

            c(l2fnorm, l2entrophy, l2normalized):=metric_cov(naive_esti, true_cov)
            output_metrics <- cbind(rbind(l2fnorm,l2entrophy,l2normalized), metric_names, d, 'naive', repeat_idx)
            error_matrix <- rbind(error_matrix, output_metrics)
            neglikeli_matrix <- rbind(neglikeli_matrix, cbind(1, truth$true_likeli, d, 'truth', repeat_idx))


            c(l2fnorm, l2entrophy, l2normalized):=metric_cov(fagqem_esti, true_cov)
            output_metrics = cbind(rbind(l2fnorm,l2entrophy,l2normalized), metric_names, d, 'fagqem', repeat_idx)
            error_matrix <- rbind(error_matrix, output_metrics)
            neglikeli_matrix <- rbind(neglikeli_matrix, cbind(1:length(efm_fagqem$like_list), efm_fagqem$like_list, d, 'fagqem', repeat_idx))
        } # end of d
    } # end of repeat_idx

    # [Create error df for analysis]
    error_matrix <- error_matrix[-1, ]; neglikeli_matrix <- neglikeli_matrix[-1,]
    error_df <- data.frame(error_matrix)
    rownames(error_df ) <- 1:nrow(error_df )
    colnames(error_df) =c('error', 'error_type', 'd', 'esti_method', 'repeat_idx')
    error_df[c('error', 'd' )] = apply(error_df[c('error', 'd' )], 2, as.numeric)


    # [Save the result]
    save_dir = paste(load_dir, 'total_466_20', '.RData', sep = '')
    save(error_matrix, file =  save_dir)
}# end of family




# [show the boxplot]
plot_error = filter(error_df, esti_method %in% c('naive', 'fagqem'), d> 16)
ggplot(plot_error) + geom_boxplot(aes(x = as.factor(d), y = error,
                                    color = esti_method)) +
  facet_wrap(~error_type, scales = "free")

# [show the avg error]
aggerror_df <- plot_error %>% group_by(d, esti_method, error_type) %>%
  summarise(avg_error = mean(error), sd_error = sd(error))

# [show the sd error]
ggplot(aggerror_df) + geom_line(aes(x = d, y = avg_error, colour = esti_method, linetype = esti_method)) +
  facet_wrap(~error_type, scales = "free")

ggplot(aggerror_df) + geom_line(aes(x = d, y = sd_error, colour = esti_method, linetype = esti_method)) +
  facet_wrap(~error_type, scales = "free")


neglikeli_df <-data.frame(neglikeli_matrix)
colnames(neglikeli_df) = c('iter', 'neglikeli', 'd', 'esti_method', 'repeat_idx')
neglikeli_df$neglikeli <- as.numeric(unlist(neglikeli_df['neglikeli']))


plot_neglikli = filter(neglikeli_df, !(esti_method %in% c('truth')))

ggplot(plot_neglikli) + geom_point(aes(x = as.numeric(iter), y = log(neglikeli),
                                       colour = as.factor(esti_method))) +
    facet_wrap(~as.factor(d), scales = "free")

