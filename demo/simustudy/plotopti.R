# EFM Optimization Results Plotting (Paper Section 4.1)
#
# This script creates plots for optimization efficiency comparison
# across different algorithms and exponential families.
#
# Prerequisites: Run ../../install_dependencies.R first
# Note: Requires results from optiexp.R runs

# Load required packages
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

# Use local results directory
results_dir <- "results"
if (!dir.exists(results_dir)) {
  stop("Results directory not found. Please run optiexp.R first to generate results.")
}



# Experiment configuration
n = 512; max_epoch = 25; d = 512; q_list <- c(6, 8, 12)
family_namelist <- c('binomial')
algo_names <- c('ps', 'sml', 'lapl', 'em')


sample_list <- c(50, 300, 500)
summary_table = data.frame( matrix(ncol = 8)); colnames(summary_table) = c('loss', 'Model', 'size', 'algo', 'time', 'd', 'q', 'comp_time')
# for (d in d_list){
for (q in q_list){
  for (family_idx in 1:length(family_namelist)){

    for (algo_idx in 1:length(algo_names)){
      if (algo_idx <=2){
        for (sample_idx in 1:length(sample_list)){
          load_name = file.path(results_dir, paste( algo_names[algo_idx],
                                    family_namelist[family_idx], paste('s', sample_list[sample_idx], sep =''),
                                    paste('d', d, sep = ''),
                                    paste('q', q, sep = ''),
                                    paste('T', max_epoch, sep= ''), sep = '_', '.RData'))
          skip_to_next <- FALSE

          tryCatch({
          load(load_name)
          opti_iter = length(efm_result$like_list)
          #efm_result$efm_time <- as.numeric(efm_result$efm_time,  units="secs")

          efm_result$efm_time <- as.double(efm_result$efm_time) + as.double(efm_result$eval_time) - as.numeric(efm_result$eval_time,  units="secs")

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
        load_name = file.path(results_dir, paste( algo_names[algo_idx],
                                  family_namelist[family_idx], paste('s', sample_list[sample_idx], sep =''),
                                  paste('d', d, sep = ''),
                                  paste('q', q, sep = ''),
                                  paste('T', max_epoch, sep= ''), sep = '_', '.RData'))
        tryCatch(
          {
            load(load_name)
            opti_iter = length(efm_result$like_list)
            #efm_result$efm_time <- as.numeric(efm_result$efm_time,  units="secs")
            efm_result$efm_time <- as.double(efm_result$efm_time) + as.double(efm_result$eval_time) - as.numeric(efm_result$eval_time,  units="secs")
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

# Filter data for analysis
summary_table = filter(summary_table,  algo %in% c('ps', 'sml', 'lapl', 'em'),
                       Model %in% c('poisson', 'binomial', 'negbinom'), size %in% c(50, 300, 500))

# Create optimization comparison plot
ggplot(summary_table) + geom_point(aes(x = as.numeric(time),
                                       y= log(as.numeric(loss)/n),
                                       shape = size,
                                       colour = algo), alpha = 0.8
                                   )+
  theme_bw() + xlab('Adam Steps') + ylab('Avged Negative likelihood') +
  facet_wrap(~Model+q, scales = "free") +
  theme(legend.position="bottom")

ggplot(summary_table) + geom_point(aes(x = as.numeric(time),
                                       y= log(comp_time),
                                       shape = size,
                                       colour = algo), alpha = 0.8
)+
  theme_bw() + xlab('Adam Steps') + ylab('Avged Negative likelihood') +
  facet_wrap(~Model + q , scales = "free") +
  theme(legend.position="bottom")

# Summary table of computation times
time_table <- filter(summary_table, Model %in% c('binomial')) %>% 
  group_by(q, size, algo) %>% 
  summarize(mean_time = mean(comp_time), .groups = 'drop')

print(time_table[order(time_table$algo),], n = 30)
