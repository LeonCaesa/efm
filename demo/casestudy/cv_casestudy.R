setwd(dirname(rstudioapi::getSourceEditorContext()$path))
if (!require("devtools")) install(devtools)
if (!require("MASS")) install(MASS)
if (!require("R.matlab")) install(R.matlab)
if (!require("tidyverse")) install(tidyverse)
if (!require("snedata")) install(snedata)
if (!require("dmf")) install_github("carvalho-research/dmf")

if (!exists("foo", mode="function")) source("util_casestudy.R")
if (!exists("foo", mode="function")) source("../../R/efm.R")


# [for ORL face]
X<- t(readMat('data/ORL_64x64.mat')$fea)
label = readMat('data/ORL_64x64.mat')$gnd

# [for Fasion Mnist] # first 60,000 instances are the training set
# fashion <- download_fashion_mnist()
# X <- as.matrix(fashion[1:2000, -c(785,786)])
# label = fashion[1:2000, 786]




n = dim(X)[1]
d = dim(X)[2]
phi_star = mean(X)^2/(sd(X)^2 - mean(X))
factor_family1 = negative.binomial(phi_star)
rank_esti = onatski_rank(X, factor_family1, q_max = d-5)

eigen_values1= eigen(cov(tcrossprod(rank_esti$L,rank_esti$V)))$value
plot_eigen = data.frame(negbinom = -diff(eigen_values1)[1:30])
ggplot(plot_eigen) + geom_point(aes(x= 1:30, y =negbinom)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Negbinom Eigen Gap') + geom_vline(xintercept = 3)+
  geom_vline(xintercept = 7)+
  theme(plot.title = element_text(hjust = 0.5))




#q = 3
#q = 7
q = 41
step_size = 0.2
#batch_size = 128
batch_size = 256
sample_size = 50

dmf_nbinom = dmf(X, factor_family1, q)
Vt_nbinom = dmf_nbinom$V; Lt_nbinom = dmf_nbinom$L
result_nbinom = batch_opti(dmf_nbinom$L, batch_size, step_size, X, factor_family = factor_family1, q,
                           max_epoch = 50,sample_size = sample_size,
                           beta1 = 0.9, beta2 = 0.999,epislon = 10^-8,
                           phi_star = phi_star, rho = 0.5, scale_weights = 1, sample_random = TRUE)


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

