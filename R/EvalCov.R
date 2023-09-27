source('./R/utils.R')
source('./R/efm.R')
library(MASS)

n_gq <- 15



names_algo <- c('lapl', 'sml', 'ps')
mse_matrix <- matrix(nrow = length(d_list), ncol =5)

#d_list = seq(10, 100, by = 10); n = 500
# factor_family = quasipoisson()
# load_dir = '/projectnb/dmfgrp/efm/CovResult/quasipoisson/'

# factor_family = negative.binomial(2)
# load_dir = '/projectnb/dmfgrp/efm/CovResult/negbinom_phi2/'

d_list = seq(16, 136, by = 20); n = 756
factor_family = negative.binomial(2)
load_dir = '/projectnb/dmfgrp/efm/CovResult/negbinom2_jianqing/'

d_count = 1
for (d in d_list){
    mse_d <- rep(0, 5)
    truth_name <- paste(paste(load_dir, 'truth', d, n, sep = '_'), '.RData', sep = '')
    load(truth_name)
    fagqem_name <- paste(paste(load_dir, 'fageqm', d, n, sep = '_'), '.RData', sep = '')
    load(fagqem_name)

    true_cov <- marg_var(truth$center, truth$V, family = factor_family, ngq = 15, dispersion = truth$phi)
    fagqem_esti <- marg_var(efm_fagqem$center, efm_identify(efm_fagqem$V), family = factor_family, ngq = 15, dispersion = efm_fagqem$dispersion)
    naive_esti <- cov(truth$X)

    mse_naive = mean((cov2cor(naive_esti) - cov2cor(true_cov)^2))
    mse_fagqem = mean( (cov2cor(fagqem_esti) - cov2cor(true_cov) )^2 )

    mse_d[1] <- mse_naive
    mse_d[2] <- mse_fagqem

    count_algo = 3
    for(algo_ in names_algo){
        efm_name <- paste(paste(load_dir, 'efm', algo_, d, n, sep = '_'), '.RData', sep = '')
        load(efm_name)

        tryCatch(
          # This is what I want to do...
          {
            efm_esti <- marg_var(efm_result$center, efm_identify(efm_result$V), family = factor_family, ngq = 15, dispersion = efm_result$dispersion)
          },
          # ... but if an error occurs, tell me what happened:
          error=function(error_message) {
            message("This is my custom message.")
          })



        mse_efm = mean((cov2cor(efm_esti) - cov2cor(true_cov) )^2)
        mse_d[count_algo] <- mse_efm
        count_algo = count_algo + 1
    }
    mse_matrix[d_count, ] <- mse_d
    d_count = d_count + 1
}

colnames(mse_matrix) <- c('naive', 'fqgqem', names_algo)

mse_matrix


# efm_esti <- marg_var(efm_result$center, efm_identify(efm_result$V), family = factor_family, ngq = 15, dispersion = efm_result$dispersion)

#

# c(mse_naive, mse_efm, mse_fagqem)
