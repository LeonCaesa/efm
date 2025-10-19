# Load the efm package (works whether installed or in development)
if (requireNamespace("efm", quietly = TRUE)) {
  library(efm)
} else {
  # For development: load from source (assuming we're in demo/casestudy/)
  devtools::load_all("../..")
}

# Load required packages
if (!require("MASS")) install.packages("MASS")
if (!require("R.matlab")) install.packages("R.matlab")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("mvtnorm")) install.packages("mvtnorm")

# Optional packages
if (!require("snedata")) {
  message("Package 'snedata' not available - some functionality may be limited")
}
if (!require("dmf")) {
  message("Package 'dmf' not available - install with: devtools::install_github('carvalho-research/dmf')")
}

# Source utility functions if they exist
if (file.exists("utilities.R")) {
  source("utilities.R")
}


# [for ORL face]
X<- readMat('data/ORL_64x64.mat')$fea
label = readMat('data/ORL_64x64.mat')$gnd
# X<- readMat('data/ORL_32x32.mat')$fea
# label = readMat('data/ORL_32x32.mat')$gnd


n = dim(X)[1]
d = dim(X)[2]
phi_star = mean(X)^2/(sd(X)^2 - mean(X))
factor_family1 = negative.binomial(phi_star)

# rank_esti = onatski_rank(X, factor_family1, q_max = d-5)
# save(rank_esti, file = 'orl_rank_1024_0408_2025_lapl_unrotated.RData')

#
# eigen_values1= eigen(cov(tcrossprod(rank_esti$L,rank_esti$V)))$value
# plot_eigen = data.frame(negbinom = -diff(eigen_values1)[1:30])
# ggplot(plot_eigen) + geom_point(aes(x= 1:30, y =negbinom)) +
#   xlab('q') +ylab('eigen diff')+
#   ggtitle('Negbinom Eigen Gap') + geom_vline(xintercept = 3)+
#   geom_vline(xintercept = 7)+
#   theme(plot.title = element_text(hjust = 0.5))




# q = 3
q = 7
# q = 41
batch_size = 256
sample_size = 300

sample_control = sample.control(sample_size = sample_size, eval_size = 500)
adam_control = adam.control(
      max_epoch = 50,
      batch_size = batch_size,
      step_size = 0.1,
      rho = 0,
      abs_tol = 1e-6,
      beta1 = 0.9,
      beta2 = 0.999,
      epsilon = 10^-8)


svd <- svd(X)
start_point = list(Vt = svd$v[,1:q], phi = 1, center = matrix(0, ncol = d))


result_nbinom = efm(X, factor_family = factor_family1, rank = q, weights = 1, start = start_point,
                    algo= 'lapl', adam_control = adam_control, sample_control = sample_control,
                    eval_likeli = TRUE)

# save(result_nbinom, file = 'orl_1024_0410_2025_rank3.RData')
# save(result_nbinom, file = 'orl_1024_0410_2025_rank7.RData')

# save(result_nbinom, file = 'orl_1024_0410_2025_rank7_unrotated_lapl.RData')
# save(result_nbinom, file = 'orl_1024_0410_2025_rank41_unrotated_lapl.RData')

# save(result_nbinom, file = 'orl_4096_0415_2025_rank41_unrotated_lapl.RData')
save(result_nbinom, file = 'orl_4096_0415_2025_rank7_unrotated_lapl.RData')
plot(result_nbinom$like_list,col = 'blue')




# [random face sampling]
plot_fit <- function(mu_hat, rnum_pixel, cnum_pixel, num_pic){
  pixels_gathered <- mu_hat %>%
    mutate(instance = row_number()) %>%
    gather(pixel, value, -instance) %>%
    tidyr::extract(pixel, "pixel", "(\\d+)", convert = TRUE) %>%
    mutate(pixel = pixel - 1,
           x = pixel %% rnum_pixel,
           y = cnum_pixel - pixel %/% rnum_pixel)
  theme_set(theme_light())

  pixels_gathered %>%
    filter(instance <= num_pic) %>%
    ggplot(aes(x, y)) + geom_raster(aes(fill= value))+
    facet_wrap(~ instance)+scale_fill_gradient(low="black",high="white")
}

rnum_pixel = 32
cnum_pixel = 32
num_pic = 1


mu_hat = tcrossprod(result_nbinom$L, result_nbinom$V)
plot_fit(as.tibble(X[1:6,]), 28, 28, 6)
plot_fit(as.tibble(factor_family1$linkinv(mu_hat)), 28, 28, 6)




# [some plotting]
#plot_negbin_emfdf = data.frame(dmf_center(result_nbinom)$L)
plot_negbin_emfdf = data.frame(result_nbinom$L)
plot_negbin_emfdf$label = as.factor(label)

library(plotly)
plot_ly(plot_negbin_emfdf, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))



# [out of sample]
fashion = X
test_data  = tail(fashion, 1000)
test_x = as.matrix(test_data)
test_label = tail(label, 1000)

#test_x = as.matrix(test_data[, -c(785,786)])
#test_label = as.matrix(test_data[, 786])

crop_x = test_x[1:10,]
set.seed(5)
#crop_x[,sample(1:784, 784/3)]=0
#crop_x[, 1:100]=0

pic_index = 5

# [compute the statistics]
L_mle = batch_mle(crop_x, result_nbinom$V, factor_family1, q)
#L_mle = batch_mle(test_x, result_nbinom$V, factor_family1, q)
L_pos = comput_mupos(L_mle[pic_index,], result_nbinom$V, factor_family1)
CholSigma = comput_CholSigma(t(L_pos), result_nbinom$V, factor_family1, scale_weights = 1)

# [simualte]
L_sim = simu_pos(L_pos, CholSigma, 100)
mu_sim = factor_family1$linkinv(tcrossprod(t(L_sim), result_nbinom$V))

library(gridExtra)
#g1 = plot_fit(as.tibble(t(test_x[pic_index,])),  19, 19, 1) + ggtitle('Original Picture') + theme(plot.title = element_text(hjust = 0.5))
g1 = plot_fit(as.tibble(t(crop_x[pic_index,])),  19, 19, 1) + ggtitle('Original Picture') + theme(plot.title = element_text(hjust = 0.5))
g2 = plot_fit(as.tibble(mu_sim),  19, 19, 6)  + ggtitle('Posterior Simulation') + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g1,g2, ncol =2)

# [------------------check quantized metrics------------------]
load("/projectnb/dmfgrp/efm/SavedExps/mnist_rank3.RData")
library(Rtsne)
library(fpc)
fit_tsne = Rtsne(X, perplexity=50, theta=0, dims=3, pca = FALSE)


num_label <- as.numeric(label)
plot_nmfdf <- data.frame(cbind(plot_nmfdf, num_label))
plot_negbin_dmfdf <- data.frame(cbind(dmf_nbinom$L, num_label))
plot_tsnedf <- data.frame(cbind(fit_tsne$Y, num_label))
plot_negbin_emfdf <- data.frame(cbind(result_nbinom$L, num_label))
colnames(plot_nmfdf)<- c('X1', 'X2', 'X3', 'label')
colnames(plot_negbin_dmfdf)<- c('X1', 'X2', 'X3', 'label')
colnames(plot_tsnedf)<- c('X1', 'X2', 'X3', 'label')
colnames(plot_negbin_emfdf)<- c('X1', 'X2', 'X3', 'label')


CH_tsne = round(calinhara(plot_tsnedf[,-4], num_label), digits = 2)
CH_nmf = round(calinhara(plot_nmfdf[,-4], num_label), digits = 2)
CH_negbin_dmf = round(calinhara(plot_negbin_dmfdf[,-4], num_label), digits = 2)
CH_negbin_efm = round(calinhara(plot_negbin_emfdf[,-4], num_label), digits = 2)


library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
library(HandTill2001)
plot_lists <- list(plot_nmfdf, plot_negbin_dmfdf, plot_tsnedf, plot_negbin_emfdf)
model_names <- c('nmf', 'negbin_dmf', 'tsne', 'negbin_efm')
for (plot_idx in 1:4){

    plot_modeldf <- plot_lists[[plot_idx]][,-4]
    splitflag = TrainTest_Flag(plot_modeldf, seed_ = 0, train_ratio = 0.5)
    plot_modeldf$y = num_label
    splitdata = list(train = plot_modeldf[splitflag$train, ], test = plot_modeldf[splitflag$test, ])


    tree_fit = Tree_tuned(splitdata$train)
    multi_fit = multinom(y~., splitdata$train, trace = FALSE)
    knn_fit = knn3(y~., splitdata$train, k = 9)

    tree_pred = predict(tree_fit, newdata = splitdata$test)
    multi_pred = predict(multi_fit, newdata = splitdata$test, type= 'prob')
    knn_pred = predict(knn_fit, newdata = splitdata$test)

    response = splitdata$test$y
    htauc_tree = auc(multcap(response = factor(response),predicted = tree_pred))
    htauc_multi = auc(multcap(response = factor(response), predicted = multi_pred))
    htauc_knn = auc(multcap(response = factor(response),predicted = knn_pred))
    ht_scores = c(htauc_tree, htauc_multi, htauc_knn)
    print(c(model_names[plot_idx], ht_scores))
}


library(plotly)
plot_ly(plot_tsnedf, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~ as.factor(label),
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
