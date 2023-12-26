source('./R/utils.R')
source('./R/efm.R')
library(MASS)
library(lemon)
library(tidyverse)


#file_name = '_total_466.RData'
file_name = 'total_466_20.RData'
file_dir = '/projectnb/dmfgrp/efm/CovResult1209/'

family_list <- c('quasipoisson', 'negbinom(20)', 'poisson', 'binomial')


agg_df = data.frame(matrix(ncol = 6))
colnames(agg_df) <- c("error", "error_type", "d" , "esti_method", "repeat_idx","family")
for (family_name in family_list){
  load_dir = paste(file_dir, family_name, '/', file_name, sep = '')
  load(load_dir)
  error_df$family = family_name

  agg_df = rbind(agg_df, error_df)

}
agg_df <- agg_df[-1, ]
agg_df$esti_method[agg_df$esti_method == 'fagqem'] <- 'efm(em)'


#plot_df <- filter(agg_df, error_type %in% c('l2normalized'), d>16)
plot_df <- filter(agg_df, error_type %in% c('l2entrophy', 'l2normalized', "l2frobenius"), d>16)
# ggplot(plot_df) +  theme_bw() +
#   geom_boxplot(aes(x = as.factor(d), y = error, color = esti_method)) +
#   ylab('error') + xlab('d') +
#   facet_wrap(~error_type + family, labeller = label_wrap_gen(multi_line=FALSE), scales = "free", nrow = 3,
#              as.table = FALSE) +
#   theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")



# png(filename = '/projectnb/dmfgrp/efm/cov_error.png', width = 10, height = 6, units = 'in', res = 300)
ggplot(plot_df) +  theme_bw() +
  geom_boxplot(aes(x = as.factor(d), y = error, color = esti_method)) +
  ylab('cov error') + xlab('d') +
  facet_rep_grid(error_type ~ family, repeat.tick.labels = TRUE, scales = 'free_y') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
# dev.off()
