library(efm)
if (!require("multinet")) install(multinet)
if (!require("tidyverse")) install(tidyverse)
if (!require("plotly")) install(plotly)
if (!require("mvtnorm")) install(mvtnorm)


#setwd('/projectnb2/dmfgrp/efm')
#if (!exists("foo", mode="function")) source('../R/utils.R')
if (!exists("foo", mode="function")) source("util_casestudy.R")

# library(mvtnorm)
# library(matrixStats)
# library(MASS)




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
