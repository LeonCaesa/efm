setwd(dirname(rstudioapi::getSourceEditorContext()$path))
if (!require("multinet")) install(multinet)
if (!require("tidyverse")) install(tidyverse)
if (!require("plotly")) install(plotly)
if (!require("mvtnorm")) install(mvtnorm)
if (!require("matrixStats")) install(matrixStats)


if (!exists("foo", mode="function")) source("util_casestudy.R")
if (!exists("foo", mode="function")) source('../../R/utils.R')
if (!exists("foo", mode="function")) source('../../R/efm.R')




net <- ml_aucs()

get_layeradj <- function(layerA, total_ids){
  n_actors = length(total_ids)
  A = matrix(0, nrow = n_actors, ncol = n_actors)
  layer_names = rownames(layerA)


  nonzeros_idx = which(as.matrix(layerA)!=0, arr.ind = TRUE)
  n_ = dim(nonzeros_idx)[1]

  for(i in 1:n_){

    from_idx = c(nonzeros_idx[i,][1])
    to_idx = c(nonzeros_idx[i,][2])

    A_fromidx = which(total_ids ==layer_names[from_idx])
    A_toidx = which(total_ids ==layer_names[to_idx])
    A[A_fromidx, A_toidx] = layerA[from_idx, to_idx]


  }
  return(A)
}


get_totaladj <- function(network, actor_list){
  A_totl = get.adjacency(as.igraph(network))
  str_order = str_sort(rownames(A_totl), numeric = TRUE)
  A_totl = A_totl[str_order, str_order]

  node_ids = rownames(A_totl)
  A_list = c()

  layer_namelist = layers_ml(network)
  n_layers = length(layer_namelist)

  for(layer_idx in 1:n_layers){
    layer_name = layer_namelist[layer_idx]
    network_layer = as.igraph(network, layers = layer_name)
    A_layer = get.adjacency(network_layer)
    A_list[[layer_idx]] = get_layeradj(A_layer, node_ids)
  }
  return(A_list)
}

get_totalweights <- function(A_list, diag0 = TRUE, multiplier = 1){
  n_actor = dim(A_list[[1]])
  n_layer = length(A_list)
  weight_matrix = matrix(0, nrow = n_actor, ncol = n_actor)
  for(i in 1:n_layer){
    A_eigen = eigen(A_list[[i]])
    weight_matrix = weight_matrix + 1/ Re(A_eigen$values[1]) * A_list[[i]]
  }
  weight_matrix[weight_matrix!=0] = weight_matrix[weight_matrix!=0]/ min(weight_matrix[weight_matrix!=0]) * multiplier
  weight_matrix[weight_matrix ==0] = min(weight_matrix[weight_matrix!=0])
  if (diag0){diag(weight_matrix) = 0}

  return(weight_matrix)
}

# # [get aggregated adj matrix]
network_totl =  as.igraph(net)
A_totl = get.adjacency(network_totl)
str_order = str_sort(rownames(A_totl), numeric = TRUE) # the nodes needs to be ordered to match with the labels
A_totl = A_totl[str_order, str_order]
A_nonzeros = A_totl[A_totl!=0]
A_totl[A_totl>1] = 1
nodes = rownames(A_totl)
A_totl = as.matrix(A_totl)




# [get label]
attributes_ml(net)
labels = unlist(get_values_ml(net, "group", actors=data.frame(actor=nodes)))
labels[labels=="G2/G3"] = "G3"
labels[labels=="G2/G6"] = "G2"
NA_Flag = labels =='NA'
n = 61


# [get agg add]
A_list = get_totaladj(net)
n = dim(A_list[[1]])[1]; p = dim(A_list[[1]])[2]
glm_weights = get_totalweights(A_list, diag = TRUE, multiplier = 3)

# [efm fit]
glm_family = binomial('logit')
rank_ = 3
PC_names = paste('PC', seq(1:rank_), sep ='')
A_dmf = as.matrix(A_totl>=1)
#A_dmf = A_totl/glm_weights

# glm_weights = matrix(rep(rowSums(A_dmf), n), nrow =n)
diag(glm_weights)= 0

adam_control = adam.control(max_epoch = 1000, batch_size = 32,
                            step_size = 0.5, rho =0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10 ^ -8)
sample_control = sample.control(sample_size = 500, eval_size = 1500)


Vstart = svd(A_totl, nu = rank_, nv = rank_)$v
center_start = rep(0, p)
dispersion_start = 1
start = list(Vt = Vstart, center = center_start, ddispersion = dispersion_start)
L_prior <- list(mean = rep(0, rank_),
                precision =  rep(1,rank_))



efm_fit = efm(A_dmf, factor_family = glm_family, rank = rank_, weights = glm_weights,
              algo = 'lapl', start = start, adam_control = adam_control, lambda_prior = L_prior,
              sample_control = sample_control, eval_likeli = TRUE)



#save(efm_fit, file = '/projectnb/dmfgrp/Laplacian_EFM/Result/CVFit/EFM_AUCS1000Epoch.RData')
#load('/projectnb/dmfgrp/Laplacian_EFM/Result/CVFit/EFM_AUCS200Epoch.RData')
#plot(efm_fit$like_list)

PCs_WDMF= efm_identify(efm_fit$V)
PCs = PCs_WDMF
# [agg and plot the results]
plot_df = data.frame(cbind(PCs, labels))[NA_Flag == FALSE,]
colnames(plot_df) = c(PC_names, 'label')
plot_df[,PC_names] = apply(plot_df[,PC_names], 2, as.numeric)



m <- list(
  l = 0,
  r = 0,
  b = 0,
  t = 50,
  pad = 4
)
plot_ly(data = plot_df, x= ~PC1, y =~ PC2, z =~PC3,  opacity= 1,
        color = ~ as.factor(label),
        type="scatter3d", mode="markers",
        marker = list(size = 5)) %>%
  layout(title = 'Weighted EFM Embedding',
         legend = list(orientation = 'h', xanchor = "center", x = 0.5), margin = m)







# [plot individual network to show sparsity]
library(fields)
library("latex2exp")
net <- ml_aucs()
network_totl <-  as.igraph(net)
A_totl <- get.adjacency(network_totl)
A_totl[A_totl!=0] <- 0



labels = get_values_ml(net, "group", actors=data.frame(actor= colnames(A_totl)))
labels[labels=="G2/G3"] = "G3"
labels[labels=="G2/G6"] = "G2"
NA_Flag = labels =='NA'

to_sort_nodes = cbind(labels$group, colnames(A_totl))
sorted_nodes = to_sort_nodes[order(to_sort_nodes[,1]),]

colnames(A_totl) <- sorted_nodes[,2]; rownames(A_totl) <-sorted_nodes[,2]
A1 <- A_totl; A2 <- A_totl; A3 <- A_totl; A4 <- A_totl; A5 <- A_totl


from_actor <- edges_ml(net)$from_actor
to_actor <- edges_ml(net)$to_actor
layer_idx <- edges_ml(net)$to_layer
for (i in 1:length(from_actor)){
  if (layer_idx[i] == "work"){
    A1[from_actor[i], to_actor[i]] <- 1
  }else if(layer_idx[i] == "coauthor"){
    A2[from_actor[i], to_actor[i]] <- 1
  }else if(layer_idx[i] == "lunch"){
    A3[from_actor[i], to_actor[i]] <- 1
  }else if(layer_idx[i] == "facebook"){
    A4[from_actor[i], to_actor[i]] <- 1
  }else if(layer_idx[i]== "leisure"){
    A5[from_actor[i], to_actor[i]] <- 1
  }
}
g1 <- graph_from_adjacency_matrix(A1, mode='undirected')
g2 <- graph_from_adjacency_matrix(A2, mode='undirected')
g3 <- graph_from_adjacency_matrix(A3, mode='undirected')
g4 <- graph_from_adjacency_matrix(A4, mode='undirected')
g5 <- graph_from_adjacency_matrix(A5, mode='undirected')


# png("/Users/caesa/Desktop/BU PhD/Dissertation/CaseStudy/figures/sparsity_network.png",
#     units="in", width=12, height=4, res=300)
par(mfrow=c(1,5), mar=c(5,1,5,1))
plot(g1, layout=layout.sphere, main="work")
plot(g2, layout=layout.sphere, main="coauthor")
plot(g3, layout=layout.sphere, main="lunch")
plot(g4, layout=layout.sphere, main="facebook")
plot(g5, layout=layout.sphere, main="leisure")
# dev.off()


image.real <- function(mat, main_name = 'NA') {
  mat <- t(mat)[,nrow(mat):1]
  #image.plot(mat, axes = FALSE, main = main_name, xaxt= "n", yaxt= "n")
  image(mat, axes = FALSE, main = main_name,  col = c("white", "black"), cex = 10)
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat),
       tick = FALSE, las = 2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat),
       tick = FALSE , las = 2)
  box()
}


#png("/projectnb/dmfgrp/efm/figures/sparsity_network2.png",
 #   units="in", width=15, height=4, res=300)
par(mfrow=c(1,5), mar=c(5,3,5,0.5))
image.real(as_adjacency_matrix(g1, sparse = FALSE), main_name = TeX('\\textbf{$A^{(1)}$-work}'))
image.real(as_adjacency_matrix(g2, sparse = FALSE), main_name = TeX('\\textbf{$A^{(2)}$-coauthor}'))
image.real(as_adjacency_matrix(g3, sparse = FALSE), main_name = TeX('\\textbf{$A^{(3)}$-lunch}'))
image.real(as_adjacency_matrix(g4, sparse = FALSE), main_name = TeX('\\textbf{$A^{(4)}$-facebook}'))
image.real(as_adjacency_matrix(g5, sparse = FALSE), main_name = TeX('$\\textbf{A^{(5)}$-leisure}'))
#dev.off()


