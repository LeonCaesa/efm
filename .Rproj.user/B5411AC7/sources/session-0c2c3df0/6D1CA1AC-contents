# [for plotting]
library(tidyverse)
library(gridExtra)
library(latex2exp)

ThreeD = FALSE
if (ThreeD){
  #load_dir = '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc_chrRecorded.RData'  
  load_dir = '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/class_result/auc_chrRecorded.RData'  
  bench_dir = "/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc_AllChr.RData"
  model_factor = 'dmf3D_phi0.1'
  comp_factor = 'PCA3D'
  
}else{
  # load_dir = '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc20_chrRecorded.RData'
  load_dir = '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/class_result/auc20_chrRecorded.RData'  
  bench_dir = "/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc20_AllChr.RData"
  model_factor = 'dmf20D_phi0.1'
  comp_factor = 'PCA20D'
}
load(load_dir)
load(bench_dir)
auc_long$model[auc_long$model == 'PCA20sD']  = 'PCA20D';auc_agg = auc_long; auc_agg$chr_num = 'all'
bench_plot = auc_agg%>% group_by(y_label, model, class_model) #%>% summarise(mean_ = mean(value))
bench_plot$chr_num = 'all'


auc_df[,1:3] = mutate_all(auc_df[,1:3], function(x) as.numeric(as.character(x)))
colnames(auc_df) = c("tree" ,"multi", "knn",   "y_label" ,  "model"  ,   "chr_num" ,  "repeats"  )
auc_df$chr_num = as.factor(as.numeric( auc_df$chr_num))
auc_long <- auc_df %>%
  gather(class_model, value, tree:knn)
auc_subdf3D = filter(auc_long, model %in% c(model_factor, comp_factor))


auc_subdf3D$y_label = recode(auc_subdf3D$y_label, super_pop = "continent", pop = "nationality")
bench_plot$y_label = recode(bench_plot$y_label, super_pop = "continent", pop = "nationality")


auc_subdf3D = rbind(bench_plot,auc_subdf3D) 

auc_subdf3D$chr_num = factor(auc_subdf3D$chr_num, levels = c('all', 1:22))

#png('/projectnb/dmfgrp/GeneDMF/figure/3DAUC.png', units="in", width=8, height=6, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/20DAUC.png', units="in", width=8, height=6, res=300)
ggplot(auc_subdf3D) + geom_boxplot( outlier.size = 0.1, aes(x = chr_num, y= value, colour = model)) +
  #geom_boxplot(data = bench_plot,  aes(x = 'Agg',y = value, colour = model), linetype = 'dotdash') +
  facet_wrap(~ class_model + y_label, ncol = 2,
             labeller = labeller(
               .multi_line = FALSE), scales = "free"
  ) + theme_bw() + theme(legend.position="bottom") + 
  ylab('Multi-AUC Score') + xlab('Chr Num') + 
  #scale_colour_manual('3D Factr Model', values=c(1,2),labels = c('DMF(Negbin)', 'PCA'))
  scale_colour_manual('20D Factr Model', values=c(1,2), labels =c('DMF(Negbin)', 'PCA'))
# dev.off()

#grid.arrange(g1,g2,nrow=2)