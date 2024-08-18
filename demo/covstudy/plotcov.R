setwd(dirname(rstudioapi::getSourceEditorContext()$path))
if (!require("tidyverse")) install(tidyverse)
if (!require("MASS")) install(MASS)
if (!require("lemon")) install(lemon)

if (!exists("foo", mode="function")) source("../../R/efm.R")
if (!exists("foo", mode="function")) source("../../R/utils.R")


#file_name = '_total_466.RData'
#file_name = 'total_466_20.RData'
#file_dir = '/projectnb/dmfgrp/efm/CovResult1209/'

family_list <- c('quasipoisson', 'negbinom(20)', 'poisson', 'binomial')


agg_df = data.frame(matrix(ncol = 6))
colnames(agg_df) <- c("error", "error_type", "d" , "esti_method", "repeat_idx","family")
for (family_name in family_list){
  load_dir = paste(file_dir, family_name, '/', file_name, sep = '')
  load(load_dir)
  error_df = data.frame(error_matrix, row.names = NULL)
  error_df$family = family_name
  colnames(error_df) <- colnames(agg_df)
  agg_df = rbind(agg_df, error_df)

}
agg_df <- agg_df[-1, ]
agg_df$esti_method[agg_df$esti_method == 'fagqem'] <- 'efm(em)'


#plot_df <- filter(agg_df, error_type %in% c('l2normalized'), d>16)
plot_df <- filter(agg_df, error_type %in% c('l2entrophy', 'l2normalized', "l2frobenius"), d>16)
plot_df$error_type[plot_df$error_type == 'l2entrophy'] = 'l2entropy'
plot_df$d =factor(plot_df$d, levels = c(66, 166, 216, 266, 316, 366, 416, 466))
# ggplot(plot_df) +  theme_bw() +
#   geom_boxplot(aes(x = as.factor(d), y = error, color = esti_method)) +
#   ylab('error') + xlab('d') +
#   facet_wrap(~error_type + family, labeller = label_wrap_gen(multi_line=FALSE), scales = "free", nrow = 3,
#              as.table = FALSE) +
#   theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")



#png(filename = '/projectnb/dmfgrp/efm/cov_error.png', width = 10, height = 6, units = 'in', res = 300)
ggplot(plot_df) +  theme_bw() +
  geom_boxplot(aes(x = as.factor(d), y = as.numeric(error), color = esti_method)) +
  ylab('cov error') + xlab('d') +
  ggh4x::facet_nested('error type' + error_type ~ 'efm family' + family, scales = 'free_y', independent = "y") +
  #facet_wrap(error_type ~ family,  scales = 'free_y') +
  #facet_rep_wrap(error_type ~ family, repeat.tick.labels = TRUE, scales = 'free') +
  theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
#dev.off()

