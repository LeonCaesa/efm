# EFM Multiplex Network Analysis (Paper Section 4.4)
#
# This script demonstrates EFM application to multiplex network data
# using the AUCS dataset for social network analysis.
#
# Prerequisites: Run ../../install_dependencies.R first

# Load EFM package
devtools::load_all("../..")

# Load required packages
if (!require("multinet")) install.packages("multinet")
if (!require("plotly")) install.packages("plotly")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("stringr")) install.packages("stringr")
if (!require("igraph")) install.packages("igraph")

# Source utility functions
if (file.exists("utilities.R")) {
  source("utilities.R")
}




# Load AUCS multiplex network dataset
net <- ml_aucs()

# Process aggregated adjacency matrix
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

cat("EFM fitting completed. Final loss:", tail(efm_fit$like_list, 1), "\n")

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




# Create results directory
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE)
}

# Create network sparsity visualization
png(file.path(results_dir, "sparsity_network.png"),
    units="in", width=15, height=4, res=300)

# Use plotmath for mathematical expressions
title_expr <- function(k, label) bquote(bold(A^{.(k)}) ~ "--" ~ .(label))
TITLE_SIZE <- 2
par(mfrow = c(1,5), mar = c(5,3,5,0.5))
image.real(as_adjacency_matrix(g1, sparse = FALSE),
           main_name = title_expr(1, "Work"),      cex_main = TITLE_SIZE)
image.real(as_adjacency_matrix(g2, sparse = FALSE),
           main_name = title_expr(2, "Coauthor"),  cex_main = TITLE_SIZE)
image.real(as_adjacency_matrix(g3, sparse = FALSE),
           main_name = title_expr(3, "Lunch"),     cex_main = TITLE_SIZE)
image.real(as_adjacency_matrix(g4, sparse = FALSE),
           main_name = title_expr(4, "Facebook"),  cex_main = TITLE_SIZE)
image.real(as_adjacency_matrix(g5, sparse = FALSE),
           main_name = title_expr(5, "Leisure"),   cex_main = TITLE_SIZE)
dev.off()

cat("Network sparsity plot saved to:", file.path(results_dir, "sparsity_network.png"), "\n")
