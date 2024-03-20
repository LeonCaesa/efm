library(tidyverse)

#file_dir = '/projectnb/dmfgrp/efm/OptiResult1221/'
#family_namelist <- c('negbinom', 'poisson', 'Gamma')
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1223/'
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1225/' #init true + rnorm
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1226/' #family_init + rnorm(1)
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0104/' #family_init + rnorm(0.3); alpha = 0.05
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0106/' #family_init + rnorm(0.5); alpha = 0.5(d=5), 0.05(d=10)
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0107/' #family_init + rnorm(0.5); alpha=0.1, Vsigma = 0.5
file_dir = '/projectnb/dmfgrp/efm/OptiResult0108/' #family_init + rnorm(0.5); alpha=0.05, Vsigma = 0.3
family_namelist <- c('negbinom', 'poisson', 'binomial')
algo_names <- c('ps', 'sml', 'em', 'lapl')
sample_list <- c(50, 300, 500)
d_list <- c(512)


n = 512; max_epoch = 25; q =2
summary_table = data.frame( matrix(ncol = 7)); colnames(summary_table) = c('loss', 'Model', 'size', 'algo', 'time', 'd', 'q')
for (d in d_list){
  for (family_idx in 1:3){

    for (algo_idx in 1:length(algo_names)){
      if (algo_idx <=2){
        for (sample_idx in 1:3){
          load_name = paste( file_dir,
                             paste( algo_names[algo_idx],
                                    family_namelist[family_idx], paste('s', sample_list[sample_idx], sep =''),
                                    paste('d', d, sep = ''),
                                    paste('q', q, sep = ''),
                                    paste('T', max_epoch, sep= ''), sep = '_'),
                             '.RData', sep ='')
          load(load_name)
          opti_iter = length(efm_result$like_list)
          efm_result$efm_time <- as.numeric(efm_result$efm_time,  units="secs")
          #time_unit = cumsum(1:opti_iter * efm_result$efm_time/ opti_iter)
          time_unit = 1:length(efm_result$like_list)
          temp_row = cbind( efm_result$like_list, family_namelist[family_idx],
                            sample_list[sample_idx], algo_names[algo_idx],
                            time_unit, d, q)
          colnames(temp_row) = colnames(summary_table)
          summary_table = rbind(summary_table, temp_row)
        } # end of sample_idx
      }else{
        sample_idx = 1
        load_name = paste( file_dir,
                           paste( algo_names[algo_idx],
                                  family_namelist[family_idx], paste('s', sample_list[sample_idx], sep =''),
                                  paste('d', d, sep = ''),
                                  paste('q', q, sep = ''),
                                  paste('T', max_epoch, sep= ''), sep = '_'),
                           '.RData', sep ='')
        tryCatch(
          {
            load(load_name)
            opti_iter = length(efm_result$like_list)
            efm_result$efm_time <- as.numeric(efm_result$efm_time,  units="secs")
            #time_unit = cumsum(1:opti_iter * efm_result$efm_time/ opti_iter)
            time_unit = 1:length(efm_result$like_list)
            temp_row = cbind( efm_result$like_list, family_namelist[family_idx],
                              sample_list[sample_idx], algo_names[algo_idx],
                              time_unit, d, q)
            colnames(temp_row) = colnames(summary_table)
            summary_table = rbind(summary_table, temp_row)
          }, error = function(e) {print(load_name)}
        )
        }# end of algo if else
      } #end of algo
    } # end of family
} # end of d


summary_table= summary_table [-1,]
summary_table$time = as.numeric(summary_table$time)


summary_table = filter(summary_table, time<=25000, algo %in% c('ps', 'sml', 'lapl', 'em'),
                       Model %in% c('poisson', 'binomial', 'negbinom'), size %in% c(50, 300, 500))
#summary_table$d = as.factor(summary_table$d)

#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep5.png', width = 8, height = 4, units = 'in', res = 300)
#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep10.png', width = 8, height = 4, units = 'in', res = 300)
png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep512.png', width = 8, height = 4, units = 'in', res = 300)
ggplot(summary_table) + geom_point(aes(x = as.numeric(time),
                                       y= log(as.numeric(loss)/n),
                                       shape = size,
                                       colour = algo), alpha = 0.8
                                   )+
  theme_bw() + xlab('Adam Steps') + ylab('Avged Negative likelihood') + facet_wrap(~Model, scales = "free") +
  theme(legend.position="bottom")
dev.off()

  #+
    # scale_shape_manual(values = c(0 ,4, 1),
    #                    labels = c('50','300','500')) +
# binomial, em/lapl needs good initialization

