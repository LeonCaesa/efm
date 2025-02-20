load("/projectnb/dmfgrp/Laplacian_EFM/Result/CVFit/ORLFace_epoch20_algo_lapl_lr0.1_b400_q40_decay_0.5.RData")

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

plot_cfit<- function(mu_hat, rnum_pixel, cnum_pixel, num_pic){
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
    pixels_gathered %>%
      filter(instance <= num_pic) %>%
      ggplot(aes(x, y)) + geom_raster(aes(fill= value))+
      facet_wrap(~ instance) + scale_fill_gradient(low="black",high="white")
  }else{
    pixels_gathered$rgb_groups = factor(pixels_gathered$rgb_groups, labels = c("R", "G", "B"))
    pixels_gathered =  pixels_gathered %>%pivot_wider(id_cols = c(instance, x, y, pixel, rgb_groups),
                                                      names_from = rgb_groups, values_from = value)%>%
      group_by(instance, x,y) %>%mutate(R = mean(R, na.rm = TRUE),
                                        G= mean(G, na.rm = TRUE),
                                        B = mean(B, na.rm = TRUE))%>%ungroup() %>%select(-pixel)%>%distinct()
    pixels_gathered %>%
      filter(instance <= num_pic) %>%
      ggplot(aes(x, y)) + geom_raster(aes(fill= rgb(R/255, G/255, B/255)))+
      scale_fill_identity()+ facet_wrap(~ instance, scales = "free")

  }
}


data_dir = '/projectnb/dmfgrp/Exponential_Factor_Model/data'
ORL_datadir = paste(data_dir, '/ORL_32x32.mat', sep ='')
X <- readMat(ORL_datadir)$fea
label = readMat(ORL_datadir)$gnd

n = dim(X)[1];p = dim(X)[2]
glm_weights = matrix(1, nrow = n, ncol = p)

rnum_pixel = 32;cnum_pixel = 32
test_idx = 2; num_pic = 4
plot_cfit(as.tibble(t(X[test_idx,])), rnum_pixel, cnum_pixel, num_pic )



L_esti <- t(mapply(ridge_coef, asplit(t(X), 1), asplit( matrix(1, nrow = p, ncol = n), 1), MoreArgs = list(Vt = efm_fit$V, factor_family = efm_fit$family)))

plot_cfit(as.tibble(t(X[test_idx,])), rnum_pixel, cnum_pixel, num_pic )

plot_cfit(as.tibble(t(L_esti[, 1])), rnum_pixel, cnum_pixel, 1)
plot_cfit(as.tibble(t(L_esti[, 2])), rnum_pixel, cnum_pixel, 1)
plot_cfit(as.tibble(t(L_esti[, 3])), rnum_pixel, cnum_pixel, 1)







