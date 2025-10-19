
# Utility functions for EFM case studies
# This file contains helper functions used by the case study scripts

# Network analysis utility functions
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

# Visualization utility functions
image.real <- function(mat, main_name = "NA", cex_main = 1.3) {
  mat <- t(mat)[, nrow(mat):1]
  image(mat, axes = FALSE, col = c("white", "black"))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat), tick = FALSE, las = 2)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat), tick = FALSE, las = 2)
  box()
  title(main = main_name, cex.main = cex_main)
}

# [classification auc score computation]
TrainTest_Flag<-function(plot_df, seed_, train_ratio =0.7){
  train_size = floor(train_ratio * nrow(plot_df))
  set.seed(seed_)
  split1<- sample(c(rep(0, train_size), rep(1, nrow(plot_df)- train_size)))
  list(train = split1==0, test = split1==1)}

TrainTest_Split<-function(plot_df, label_name, splitflag, label_level){
  if (is.element( 'PC20', colnames(plot_df))){
    # subset_cols = c(c('PC1', 'PC2', 'PC3'), label_name)}
    subset_cols = c(paste('PC', 1:20, sep = ''), label_name)}
  else{
    subset_cols = c(c('PC1', 'PC2', 'PC3'), label_name)}
  reg_df = plot_df[subset_cols]
  p = dim(reg_df)[2]; colnames(reg_df)[p] = 'y'
  reg_df$y = factor(reg_df$y, levels= label_level)

  # return (list(train = reg_df[splitflag$train, ], test = reg_df[!splitflag$test, ]))}
  return (list(train = reg_df[splitflag$train, ], test = reg_df[splitflag$test, ]))}

Tree_tuned<-function(data){
  tree_fit = rpart(y ~ ., data = data, parms = list(split ="information"),
                   control = c(cp = 0, xval =10), method= 'class')

  best_cp = tree_fit$cptable[which.min(tree_fit$cptable[, "xerror"]),"CP"]

  tree_fit = prune(tree_fit, best_cp)
  return(tree_fit)}
