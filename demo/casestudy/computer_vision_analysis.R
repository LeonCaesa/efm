cat("=== EFM Computer Vision Analysis ===\n")
cat("Starting analysis at:", format(Sys.time()), "\n\n")

# Load the efm package (works whether installed or in development)
cat("Loading EFM package...\n")
if (requireNamespace("efm", quietly = TRUE)) {
  library(efm)
} else {
  # For development: load from source (assuming we're in demo/casestudy/)
  devtools::load_all("../..")
}
cat("✓ EFM package loaded successfully\n\n")

# Load required packages
cat("Loading required packages...\n")
if (!require("glmnet")) install.packages("glmnet")
if (!require("MASS")) install.packages("MASS")
if (!require("R.matlab")) install.packages("R.matlab")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("mvtnorm")) install.packages("mvtnorm")
cat("✓ Core packages loaded\n\n")

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
cat("Loading ORL face dataset...\n")
X<- readMat('data/ORL_64x64.mat')$fea
label = readMat('data/ORL_64x64.mat')$gnd
# X<- readMat('data/ORL_32x32.mat')$fea
# label = readMat('data/ORL_32x32.mat')$gnd

cat("✓ Dataset loaded - Dimensions:", dim(X)[1], "samples x", dim(X)[2], "features\n")


n = dim(X)[1]
d = dim(X)[2]
phi_star = mean(X)^2/(sd(X)^2 - mean(X))
factor_family1 = negative.binomial(phi_star)

cat("Data preprocessing completed:\n")
cat("  - Number of samples (n):", n, "\n")
cat("  - Number of features (d):", d, "\n")
cat("  - Phi parameter:", round(phi_star, 4), "\n\n")

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

cat("Setting up EFM parameters:\n")
cat("  - Rank (q):", q, "\n")
cat("  - Batch size:", batch_size, "\n")
cat("  - Sample size:", sample_size, "\n")
cat("  - Max epochs: 50\n\n")

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

cat("Computing SVD for initialization...\n")
svd <- svd(X)
start_point = list(Vt = svd$v[,1:q], phi = 1, center = matrix(0, ncol = d))
cat("✓ Initialization completed\n\n")


cat("Starting EFM fitting...\n")
cat("This may take several minutes. Progress will be shown below:\n")
cat("Algorithm: Laplacian approximation\n")
cat("Expected runtime: 5-15 minutes depending on system\n\n")

start_time <- Sys.time()
result_nbinom = efm(X, factor_family = factor_family1, rank = q, weights = 1, start = start_point,
                    algo= 'lapl', adam_control = adam_control, sample_control = sample_control,
                    eval_likeli = TRUE)
end_time <- Sys.time()

cat("\n✓ EFM fitting completed!\n")
cat("Runtime:", round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "minutes\n")
cat("Final likelihood:", tail(result_nbinom$like_list, 1), "\n\n")

# Save results
cat("Saving results...\n")
save(result_nbinom, file = 'orl_4096_0415_2025_rank7_unrotated_lapl.RData')
cat("✓ Results saved to: orl_4096_0415_2025_rank7_unrotated_lapl.RData\n\n")

cat("Plotting likelihood convergence...\n")
plot(result_nbinom$like_list,col = 'blue')




cat("Creating visualization functions...\n")

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


cat("Generating face reconstructions...\n")
mu_hat = tcrossprod(result_nbinom$L, result_nbinom$V)
cat("Creating original faces plot...\n")
plot_fit(as.tibble(X[1:6,]), 28, 28, 6)
cat("Creating reconstructed faces plot...\n")
plot_fit(as.tibble(factor_family1$linkinv(mu_hat)), 28, 28, 6)
cat("✓ Face visualization completed\n\n")




# [some plotting]
cat("Creating 3D embedding visualization...\n")
#plot_negbin_emfdf = data.frame(dmf_center(result_nbinom)$L)
plot_negbin_emfdf = data.frame(result_nbinom$L)
plot_negbin_emfdf$label = as.factor(label)

library(plotly)
p <- plot_ly(plot_negbin_emfdf, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
print(p)
cat("✓ 3D embedding plot created\n\n")



# [out of sample]
cat("Starting out-of-sample analysis...\n")
fashion = X
test_data  = tail(fashion, 1000)
test_x = as.matrix(test_data)
test_label = tail(label, 1000)
cat("✓ Test data prepared:", nrow(test_data), "samples\n")

#test_x = as.matrix(test_data[, -c(785,786)])
#test_label = as.matrix(test_data[, 786])

crop_x = test_x[1:10,]
set.seed(5)
#crop_x[,sample(1:784, 784/3)]=0
#crop_x[, 1:100]=0

pic_index = 5

# [compute the statistics]
cat("Computing posterior statistics for image", pic_index, "...\n")
L_mle = batch_mle(crop_x, result_nbinom$V, factor_family1, q)
#L_mle = batch_mle(test_x, result_nbinom$V, factor_family1, q)
L_pos = comput_mupos(L_mle[pic_index,], result_nbinom$V, factor_family1)
CholSigma = comput_CholSigma(t(L_pos), result_nbinom$V, factor_family1, scale_weights = 1)
cat("✓ Posterior statistics computed\n")

# [simualte]
cat("Running posterior simulation...\n")
L_sim = simu_pos(L_pos, CholSigma, 100)
mu_sim = factor_family1$linkinv(tcrossprod(t(L_sim), result_nbinom$V))

library(gridExtra)
cat("Creating comparison plots...\n")
#g1 = plot_fit(as.tibble(t(test_x[pic_index,])),  19, 19, 1) + ggtitle('Original Picture') + theme(plot.title = element_text(hjust = 0.5))
g1 = plot_fit(as.tibble(t(crop_x[pic_index,])),  19, 19, 1) + ggtitle('Original Picture') + theme(plot.title = element_text(hjust = 0.5))
g2 = plot_fit(as.tibble(mu_sim),  19, 19, 6)  + ggtitle('Posterior Simulation') + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g1,g2, ncol =2)
cat("✓ Posterior simulation completed\n\n")

# [------------------check quantized metrics------------------]
cat("Starting quantitative evaluation...\n")
# Skip loading external file that may not exist
# load("/projectnb/dmfgrp/efm/SavedExps/mnist_rank3.RData")
library(Rtsne)
library(fpc)
cat("Running t-SNE for comparison...\n")
fit_tsne = Rtsne(X, perplexity=50, theta=0, dims=3, pca = FALSE)
cat("✓ t-SNE completed\n")


cat("Computing clustering metrics...\n")
num_label <- as.numeric(label)
# Skip methods that require external data
# plot_nmfdf <- data.frame(cbind(plot_nmfdf, num_label))
# plot_negbin_dmfdf <- data.frame(cbind(dmf_nbinom$L, num_label))
plot_tsnedf <- data.frame(cbind(fit_tsne$Y, num_label))
plot_negbin_emfdf <- data.frame(cbind(result_nbinom$L, num_label))
# colnames(plot_nmfdf)<- c('X1', 'X2', 'X3', 'label')
# colnames(plot_negbin_dmfdf)<- c('X1', 'X2', 'X3', 'label')
colnames(plot_tsnedf)<- c('X1', 'X2', 'X3', 'label')
colnames(plot_negbin_emfdf)<- c('X1', 'X2', 'X3', 'label')

CH_tsne = round(calinhara(plot_tsnedf[,-4], num_label), digits = 2)
# CH_nmf = round(calinhara(plot_nmfdf[,-4], num_label), digits = 2)
# CH_negbin_dmf = round(calinhara(plot_negbin_dmfdf[,-4], num_label), digits = 2)
CH_negbin_efm = round(calinhara(plot_negbin_emfdf[,-4], num_label), digits = 2)

cat("Calinski-Harabasz Index:\n")
cat("  - t-SNE:", CH_tsne, "\n")
cat("  - EFM:", CH_negbin_efm, "\n\n")


cat("Running classification evaluation...\n")
library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
library(HandTill2001)

# Only evaluate available methods
plot_lists <- list(plot_tsnedf, plot_negbin_emfdf)
model_names <- c('tsne', 'negbin_efm')

for (plot_idx in 1:2){
    cat("Evaluating", model_names[plot_idx], "...\n")
    
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
    cat("  Results:", model_names[plot_idx], "- Tree:", round(htauc_tree, 3), 
        "Multinomial:", round(htauc_multi, 3), "KNN:", round(htauc_knn, 3), "\n")
}
cat("✓ Classification evaluation completed\n\n")


cat("Creating final t-SNE visualization...\n")
library(plotly)
p_tsne <- plot_ly(plot_tsnedf, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~ as.factor(label),
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
print(p_tsne)

cat("\n=== Analysis Complete ===\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2), "minutes\n")
cat("Results saved and visualizations created successfully!\n")
