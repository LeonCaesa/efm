if (!require("tidyverse")) install(tidyverse)

#file_dir = '/projectnb/dmfgrp/efm/OptiResult1221/'
#family_namelist <- c('negbinom', 'poisson', 'Gamma')
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1223/'
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1225/' #init true + rnorm
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1226/' #family_init + rnorm(1)
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0104/' #family_init + rnorm(0.3); alpha = 0.05
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0106/' #family_init + rnorm(0.5); alpha = 0.5(d=5), 0.05(d=10)
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0107/' #family_init + rnorm(0.5); alpha=0.1, Vsigma = 0.5
#file_dir = '/projectnb/dmfgrp/efm/OptiResult0108/' #family_init + rnorm(0.5); alpha=0.05, Vsigma = 0.3


#file_dir = '/projectnb/dmfgrp/efm/OptiResult1215_2024/InitSD01/' #family_init + rnorm(0.1); alpha=0.05, Vsigma = 0.3
#file_dir = '/projectnb/dmfgrp/efm/OptiResult1215_2024/' #family_init + rnorm(0.5); alpha=0.05, Vsigma = 0.3
file_dir = '/projectnb/dmfgrp/efm/OptiResult0118_2025/' #family_init + rnorm(0.5); alpha=0.05, Vsigma = 0.3



# d_list <- c(5, 10, 512)
# n = 512; max_epoch = 25; q =2
# family_namelist <- c('negbinom', 'poisson', 'binomial')
# algo_names <- c('ps', 'sml', 'em', 'lapl')

#n = 512; max_epoch = 25; d = 512; q_list <- c(50, 100, 150, 200, 250, 300)
#n = 512; max_epoch = 25; d = 512; q_list <- c(10, 20, 30, 40, 50)
n = 512; max_epoch = 25; d = 512; q_list <- c(6, 8, 12) # choose 3
family_namelist <- c('binomial')
# #algo_names <- c('ps', 'sml', 'lapl', 'em')
algo_names <- c('ps', 'sml', 'lapl')


sample_list <- c(50, 300, 500)
summary_table = data.frame( matrix(ncol = 8)); colnames(summary_table) = c('loss', 'Model', 'size', 'algo', 'time', 'd', 'q', 'comp_time')
# for (d in d_list){
for (q in q_list){
  for (family_idx in 1:length(family_namelist)){

    for (algo_idx in 1:length(algo_names)){
      if (algo_idx <=2){
        for (sample_idx in 1:length(sample_list)){
          load_name = paste( file_dir,
                             paste( algo_names[algo_idx],
                                    family_namelist[family_idx], paste('s', sample_list[sample_idx], sep =''),
                                    paste('d', d, sep = ''),
                                    paste('q', q, sep = ''),
                                    paste('T', max_epoch, sep= ''), sep = '_'),
                             '.RData', sep ='')
          skip_to_next <- FALSE

          tryCatch({
          load(load_name)
          opti_iter = length(efm_result$like_list)
          efm_result$efm_time <- as.numeric(efm_result$efm_time,  units="secs")
          #time_unit = cumsum(1:opti_iter * efm_result$efm_time/ opti_iter)
          time_unit = 1:length(efm_result$like_list)
          temp_row = cbind( efm_result$like_list, family_namelist[family_idx],
                            sample_list[sample_idx], algo_names[algo_idx],
                            time_unit, d, q, efm_result$efm_time/length(efm_result$like_list))
          colnames(temp_row) = colnames(summary_table)
          summary_table = rbind(summary_table, temp_row)}, error = function(e) { skip_to_next <<- TRUE})
          if(skip_to_next) { next }

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
                              time_unit, d, q, efm_result$efm_time/length(efm_result$like_list))
            colnames(temp_row) = colnames(summary_table)
            summary_table = rbind(summary_table, temp_row)
          }, error = function(e) {print(load_name)}
        )
        }# end of algo if else
      } #end of algo
    } # end of family
} # end of d


summary_table= summary_table [-1,]

summary_table$comp_time = as.numeric(summary_table$comp_time)
summary_table$time = as.numeric(summary_table$time)
summary_table$q = paste('rank = ', summary_table$q, sep = '')
#summary_table$q = factor(summary_table$q,levels=c("rank = 4", "rank = 6", "rank = 8", "rank = 10", "rank = 12"))


# summary_table = filter(summary_table, time<=25000, algo %in% c('ps', 'sml', 'lapl', 'em'),
#                        Model %in% c('poisson', 'binomial', 'negbinom'), size %in% c(50, 300, 500))
# summary_table = filter(summary_table, time<=2500000, algo %in% c('ps', 'sml', 'lapl', 'em'),
#                        Model %in% c('poisson', 'binomial', 'negbinom'), size %in% c(50, 300, 500))
summary_table = filter(summary_table,  algo %in% c('ps', 'sml', 'lapl', 'em'),
                       Model %in% c('poisson', 'binomial', 'negbinom'), size %in% c(50, 300, 500))
#summary_table$d = as.factor(summary_table$d)

#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep5.png', width = 8, height = 4, units = 'in', res = 300)
#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep10.png', width = 8, height = 4, units = 'in', res = 300)
#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep512.png', width = 8, height = 4, units = 'in', res = 300)
#png(filename = '/projectnb/dmfgrp/efm/figures/EFMOptiComparep512Largeq.png', width = 8, height = 4, units = 'in', res = 300)
ggplot(summary_table) + geom_point(aes(x = as.numeric(time),
                                       y= log(as.numeric(loss)/n),
                                       shape = size,
                                       colour = algo), alpha = 0.8
                                   )+
  theme_bw() + xlab('Adam Steps') + ylab('Avged Negative likelihood') +
  facet_wrap(~Model+q, scales = "free") +
  theme(legend.position="bottom")
#dev.off()

ggplot(summary_table) + geom_point(aes(x = as.numeric(time),
                                       y= log(comp_time),
                                       shape = size,
                                       colour = algo), alpha = 0.8
)+
  theme_bw() + xlab('Adam Steps') + ylab('Avged Negative likelihood') +
  facet_wrap(~Model + q , scales = "free") +
  theme(legend.position="bottom")

time_table <- filter(summary_table, Model %in% c('binomial')) %>% group_by(q, size, algo) %>% summarize(mean_time = mean(comp_time))

print(time_table[order(time_table$algo),], n = 30)

  #+
    # scale_shape_manual(values = c(0 ,4, 1),
    #                    labels = c('50','300','500')) +
# binomial, em/lapl needs good initialization

