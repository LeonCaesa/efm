# load("/projectnb/dmfgrp/Laplacian_EFM/Result/CVFit/ORLFace_epoch20_algo_lapl_lr0.1_b400_q40_decay_0.5.RData")

# load("/projectnb/dmfgrp/efm/demo/casestudy/orl_1024_0410_2025_rank7.RData")
load('orl_1024_0410_2025_rank7_unrotated_lapl.RData')
efm_fit <- result_nbinom

rnum_pixel = 32;cnum_pixel = 32



if (!require("R.matlab")) install.packages("R.matlab")
if (!require("tidyverse")) install(tidyverse)
if (!require("glmnet")) install(glmnet)

ridge_coef <- function(X_vec, weight_vec, Vt, factor_family){
  d = dim(Vt)[1]
  sd_scalar = sqrt(var(X_vec)*(d-1)/d)
  pen_result <- glmnet(x = Vt, y= X_vec, family = factor_family, alpha = 0, lambda=1, weights = weight_vec,
                       intercept = FALSE, standardize = FALSE, thresh= 1e-10,
                       type.logistic = c("Newton"))
  as.vector(coef(pen_result, s = sd_scalar * 1/d, exact = TRUE, x = Vt, y = X_vec,
                 family = factor_family,
                 weights = weight_vec))[-1]
}

plot_cfit<- function(mu_hat, rnum_pixel, cnum_pixel, num_pic, col_row = NULL){
  pixels_gathered = mu_hat %>%
    mutate(instance = row_number()) %>%
    gather(pixel, value, -instance) %>%
    tidyr::extract(pixel, "pixel", "(\\d+)", convert = TRUE) %>%
    mutate(pixel = pixel - 1,
           x = pixel %% (rnum_pixel*cnum_pixel)%% rnum_pixel,
           y = cnum_pixel - pixel %% (rnum_pixel*cnum_pixel)%/% rnum_pixel,
           rgb_groups = factor(pixel %/% (rnum_pixel * rnum_pixel)))
  pixels_gathered = pixels_gathered%>%group_by(instance)%>%mutate(value = (value-mean(value))/sd(value))
  if (length(unique(pixels_gathered$rgb_groups)) ==1){
    p_plot <- pixels_gathered %>%
      filter(instance <= num_pic) %>%
      ggplot(aes(x, y)) + geom_raster(aes(fill= value))+
      facet_wrap(~ instance) + scale_fill_gradient(low="black",high="white")

    if (is.null(col_row)){
      p_plot <- p_plot + facet_wrap(~ instance, scales = "free")
    }else{
      p_plot <- p_plot + facet_wrap(~ instance, scales = "free",
                                    nrow = col_row[1], ncol = col_row[2])
    }

  }else{
    pixels_gathered$rgb_groups = factor(pixels_gathered$rgb_groups, labels = c("R", "G", "B"))
    pixels_gathered =  pixels_gathered %>%pivot_wider(id_cols = c(instance, x, y, pixel, rgb_groups),
                                                      names_from = rgb_groups, values_from = value)%>%
      group_by(instance, x,y) %>%mutate(R = mean(R, na.rm = TRUE),
                                        G= mean(G, na.rm = TRUE),
                                        B = mean(B, na.rm = TRUE))%>%ungroup() %>%select(-pixel)%>%distinct()
    p_plot <- pixels_gathered %>%
      filter(instance <= num_pic) %>%
      ggplot(aes(x, y)) + geom_raster(aes(fill= rgb(R/255, G/255, B/255)))+
      scale_fill_identity()

    if (is.null(col_row)){
      p_plot <- p_plot + facet_wrap(~ instance, scales = "free")
    }else{
      p_plot <- p_plot + facet_wrap(~ instance, scales = "free",
                                    nrow = col_row[1], ncol = col_row[2])
      }

  }
  return (p_plot + coord_flip() + scale_x_reverse())
}


data_dir = '/projectnb/dmfgrp/Exponential_Factor_Model/data'
ORL_datadir = paste(data_dir, '/ORL_32x32.mat', sep ='')
X <- readMat(ORL_datadir)$fea
label = readMat(ORL_datadir)$gnd

n = dim(X)[1];p = dim(X)[2]
glm_weights = matrix(1, nrow = n, ncol = p)


test_idx = 46; num_pic = 4
plot_cfit(as.tibble(t(X[test_idx,])), rnum_pixel, cnum_pixel, num_pic )



L_esti <- t(mapply(ridge_coef, asplit(t(X), 1), asplit( matrix(1, nrow = p, ncol = n), 1), MoreArgs = list(Vt = efm_fit$V, factor_family = efm_fit$family)))


#load("/projectnb/dmfgrp/efm/SavedExps/orl_face.RData")

plot_cfit(as.tibble(t(X[test_idx,])), rnum_pixel, cnum_pixel, num_pic )


plot_cfit(as.tibble(t(L_esti[, 1:10])), rnum_pixel, cnum_pixel, 1)
plot_cfit(as.tibble(t(L_esti[, 2])), rnum_pixel, cnum_pixel, 1)
plot_cfit(as.tibble(t(L_esti[, 3])), rnum_pixel, cnum_pixel, 1)
plot_cfit(as.tibble(t(L_esti[, 4])), rnum_pixel, cnum_pixel, 1)

plot_cfit(as.tibble(t(efm_fit$V[,1])), rnum_pixel, cnum_pixel, num_pic )



# [dmf]
library(dmf)
dmf_result <- dmf(t(X), family = efm_fit$family, rank = 40)
dmf_result
# [eigenface]

svd <- svd(scale(X, center= TRUE, scale = FALSE))
eigVec <- svd$v
eigVal <- svd$d/(ncol(X)-1)
#eigeVal_efm <- sort(apply(efm_fit$V, 2, norm, '2'), decreasing = TRUE)
eigeVal_dmf <- sort(apply(dmf_center(dmf_result)$L, 2, norm, '2'), decreasing = TRUE)


upto_rank <- 40
plot(eigVal[1:upto_rank]/ sum(eigVal[1:upto_rank]), ylim = c(0, 0.20))
points(eigeVal_efm[1:upto_rank]/ sum(eigeVal_efm[1:upto_rank]), col = 'red')
#points(eigeVal_dmf[1:upto_rank]/ sum(eigeVal_dmf[1:upto_rank]), col = 'red')


upto_rank <- 40
plot(cumsum(eigVal[1:upto_rank])/ sum(eigVal[1:upto_rank]))
points(cumsum(eigeVal_efm[1:upto_rank])/ sum(eigeVal_efm[1:upto_rank]), col = 'red')


# png(filename = '/projectnb/dmfgrp/efm/figures/eigen_orl.png', width = 12, height = 8, units = 'in', res = 300)
g1 = plot_cfit(as.tibble(t(eigVec[, 1:20])), rnum_pixel, cnum_pixel, 12, col_row = c(4, 3))
# dev.off()
#
# png(filename = '/projectnb/dmfgrp/efm/figures/efmeigen_orl.png', width = 12, height = 8, units = 'in', res = 300)
g2 = plot_cfit(as.tibble(t(L_esti[, 1:20])), rnum_pixel, cnum_pixel, 12, col_row = c(4, 3))
# dev.off()

#g3 = plot_cfit(as.tibble(t(result_nbinom$family$linkinv(result_nbinom$V[, 1:40]))), rnum_pixel, cnum_pixel, 40, col_row = c(8, 5))
# eta_esti = sweep(result_nbinom$V[, 1:10], 1, result_nbinom$center, '+')
eta_esti = sweep(result_nbinom$V[, 1:10], 1, 0, '+')

g3 = plot_cfit(as.tibble(t(eta_esti)), rnum_pixel, cnum_pixel, 10, col_row = c(8, 5))

library("gridExtra")
png(filename = '/projectnb/dmfgrp/efm/figures/eigen_orl.png', width = 12, height = 12, units = 'in', res = 300)
grid.arrange(g2 + ggtitle('EFMFace') + theme(plot.title = element_text(size = 15,hjust = 0.5),legend.position="none"),
            g1 + ggtitle('EigenFace') + theme(plot.title = element_text( size = 15, hjust = 0.5),legend.position="none"),
             ncol = 2)
dev.off()








