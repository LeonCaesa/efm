# EFM Multiplex Network Analysis (Paper Section 4.4) - Simplified Version
#
# This script demonstrates EFM application to multiplex network data
# using a simplified approach that avoids igraph compatibility issues.
#
# Prerequisites: Run install_packages.R first

# Load EFM package
devtools::load_all("../..")

# Load required packages
if (!require("plotly")) install.packages("plotly")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")

cat("Loading EFM multiplex network analysis (simplified version)...\n")

# Create synthetic multiplex network data for demonstration
# This simulates the AUCS dataset structure
set.seed(42)
n_nodes <- 61
n_layers <- 5

# Generate synthetic adjacency matrices for different layers
generate_layer_network <- function(n, density = 0.1) {
  A <- matrix(0, n, n)
  n_edges <- floor(n * (n-1) * density / 2)
  edges <- sample(which(upper.tri(A)), n_edges)
  A[edges] <- 1
  A <- A + t(A)  # Make symmetric
  return(A)
}

# Create layer-specific networks
A_work <- generate_layer_network(n_nodes, 0.08)
A_coauthor <- generate_layer_network(n_nodes, 0.06)
A_lunch <- generate_layer_network(n_nodes, 0.12)
A_facebook <- generate_layer_network(n_nodes, 0.10)
A_leisure <- generate_layer_network(n_nodes, 0.05)

A_list <- list(A_work, A_coauthor, A_lunch, A_facebook, A_leisure)
layer_names <- c("Work", "Coauthor", "Lunch", "Facebook", "Leisure")

# Create aggregated network
A_totl <- Reduce("+", A_list)
A_totl[A_totl > 1] <- 1  # Binary adjacency

# Generate synthetic node labels (groups)
node_labels <- sample(c("G1", "G2", "G3", "G4"), n_nodes, replace = TRUE, 
                     prob = c(0.3, 0.25, 0.25, 0.2))

cat("Synthetic network created with", n_nodes, "nodes and", n_layers, "layers\n")

# Compute weights for EFM
get_totalweights_simple <- function(A_list, diag0 = TRUE, multiplier = 3) {
  n_actor <- nrow(A_list[[1]])
  n_layer <- length(A_list)
  weight_matrix <- matrix(0, nrow = n_actor, ncol = n_actor)
  
  for(i in 1:n_layer) {
    # Simple weighting based on degree
    degrees <- rowSums(A_list[[i]])
    max_degree <- max(degrees)
    if (max_degree > 0) {
      layer_weight <- A_list[[i]] * (max_degree / (degrees + 1))
      weight_matrix <- weight_matrix + layer_weight
    }
  }
  
  weight_matrix[weight_matrix != 0] <- weight_matrix[weight_matrix != 0] / 
                                      min(weight_matrix[weight_matrix != 0]) * multiplier
  weight_matrix[weight_matrix == 0] <- min(weight_matrix[weight_matrix != 0])
  
  if (diag0) { diag(weight_matrix) <- 0 }
  return(weight_matrix)
}

# Calculate weights
glm_weights <- get_totalweights_simple(A_list, diag0 = TRUE, multiplier = 3)

# EFM fitting
glm_family <- binomial('logit')
rank_ <- 3
PC_names <- paste('PC', seq(1:rank_), sep = '')
A_dmf <- as.matrix(A_totl >= 1)

# Set up EFM parameters
adam_control <- adam.control(max_epoch = 100, batch_size = 32,
                            step_size = 0.5, rho = 0, abs_tol = 1e-6,
                            beta1 = 0.9, beta2 = 0.999, epsilon = 10^-8)
sample_control <- sample.control(sample_size = 200, eval_size = 500)

# Initialize starting values
Vstart <- svd(A_totl, nu = rank_, nv = rank_)$v
center_start <- rep(0, ncol(A_totl))
dispersion_start <- 1
start <- list(Vt = Vstart, center = center_start, ddispersion = dispersion_start)
L_prior <- list(mean = rep(0, rank_), precision = rep(1, rank_))

cat("Starting EFM fitting...\n")

# Fit EFM model
efm_fit <- efm(A_dmf, factor_family = glm_family, rank = rank_, weights = glm_weights,
               algo = 'lapl', start = start, adam_control = adam_control, 
               lambda_prior = L_prior, sample_control = sample_control, eval_likeli = TRUE)

cat("EFM fitting completed. Final loss:", tail(efm_fit$like_list, 1), "\n")

# Extract and identify principal components
PCs_WDMF <- efm_identify(efm_fit$V)
PCs <- PCs_WDMF

# Create plot data
plot_df <- data.frame(cbind(PCs, node_labels))
colnames(plot_df) <- c(PC_names, 'label')
plot_df[, PC_names] <- apply(plot_df[, PC_names], 2, as.numeric)

# Create 3D visualization
cat("Creating 3D visualization...\n")

m <- list(l = 0, r = 0, b = 0, t = 50, pad = 4)

p <- plot_ly(data = plot_df, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 1,
             color = ~as.factor(label), type = "scatter3d", mode = "markers",
             marker = list(size = 5)) %>%
  layout(title = 'EFM Multiplex Network Embedding',
         legend = list(orientation = 'h', xanchor = "center", x = 0.5), 
         margin = m)

# Create results directory
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE)
}

# Display the 3D plot (will show in RStudio or save as needed)
print(p)
cat("3D visualization created successfully\n")

# Create network sparsity visualization
png(file.path(results_dir, "network_layers_sparsity.png"),
    units = "in", width = 15, height = 4, res = 300)

par(mfrow = c(1, 5), mar = c(5, 3, 5, 0.5))

for (i in 1:length(A_list)) {
  image(A_list[[i]], axes = FALSE, col = c("white", "black"),
        main = paste("Layer", i, "-", layer_names[i]), cex.main = 1.2)
  box()
}

dev.off()
cat("Network sparsity plot saved to:", file.path(results_dir, "network_layers_sparsity.png"), "\n")

# Print summary
cat("\n=== Analysis Summary ===\n")
cat("Network size:", n_nodes, "nodes\n")
cat("Number of layers:", n_layers, "\n")
cat("EFM rank:", rank_, "\n")
cat("Final likelihood:", tail(efm_fit$like_list, 1), "\n")
cat("Node groups:", length(unique(node_labels)), "groups\n")

cat("\nAnalysis completed successfully!\n")
cat("Results saved in the 'results' directory.\n")
