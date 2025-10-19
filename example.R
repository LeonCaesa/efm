# EFM Package Example Script
# 
# This script demonstrates the EFM package functionality across multiple
# exponential families and optimization algorithms, showing the loss progression
# during optimization iterations.
# 
# Prerequisites: Run install_dependencies.R first to install required packages.

# Load the efm package from source
devtools::load_all(".")

# Load required dependencies
library(matrixStats)
library(MASS)
library(gaussquadr)

# [Comprehensive demonstration: Multiple families and algorithms]
set.seed(1)
d = 10
n = 50
q = 3

# Define test scenarios
test_scenarios <- list(
  list(
    name = "Binomial",
    dispersion_star = 1,
    factor_family = binomial(),
    factor_weights = 10
  ),
  list(
    name = "Negative Binomial", 
    dispersion_star = 20,
    factor_family = negative.binomial(20),
    factor_weights = 1
  ),
  list(
    name = "Poisson",
    dispersion_star = 1,
    factor_family = poisson(),
    factor_weights = 1
  ),
  list(
    name = "Quasi-Poisson",
    dispersion_star = 2,
    factor_family = quasipoisson(),
    factor_weights = 1
  )
)

# Algorithms to test
algorithms <- c("em", "lapl", "sml", "ps")

# Common parameters
L_prior <- list(mean = c(0, 0, 0), precision = 1/ c(1, 1, 1))
V_prior <- list(mean = c(0, 0, 0), sigma = c(0.5, 0.5, 0.5))
center_star <- rep(0, d)

cat("=== EFM PACKAGE DEMONSTRATION ===\n")
cat("Testing", length(test_scenarios), "families with", length(algorithms), "algorithms each\n\n")

for (scenario_idx in seq_along(test_scenarios)) {
  scenario <- test_scenarios[[scenario_idx]]
  
  cat("==================================================\n")
  cat("FAMILY", scenario_idx, ":", scenario$name, "\n")
  cat("==================================================\n")
  
  # Generate data for this scenario
  set.seed(1)  # Same seed for fair comparison
  truth <- generate_cov(n, d, L_prior, V_prior, center_star, 
                        family = scenario$factor_family,
                        weights = scenario$factor_weights, 
                        phi = scenario$dispersion_star)
  
  cat("Data range:", range(truth$X), "\n\n")
  
  # Print TRUE VALUES
  cat("TRUE PARAMETERS:\n")
  cat("  Center (first 3)    :", round(head(c(truth$center), 3), 4), "\n")
  cat("  Dispersion (first 3):", round(head(c(truth$phi), 3), 4), "\n")
  cat("  V0 matrix (first 3 rows):\n")
  print(round(truth$V0[1:3, ], 4))
  cat("\n")
  
  # Initialize
  Vstart <- svd(truth$X, nu = q, nv = q)$v
  center_start <- truth$center + rnorm(d, 0, 1)/10
  dispersion_start <- rnorm(d, 0, scenario$dispersion_star/3)
  dispersion_start <- dispersion_start - min(dispersion_start) + 0.1
  init <- list(Vt = Vstart, center = center_start, dispersion = dispersion_start)
  
  # Test each algorithm
  for (algo in algorithms) {
    cat("--- Algorithm:", toupper(algo), "---\n")
    
    start_time <- Sys.time()
    
    # Try to run the algorithm
    res <- tryCatch({
      if (algo == "em") {
        control <- list(maxit = 15, epsilon = 1e-6, trace = TRUE)
        efm(truth$X/scenario$factor_weights, 
            factor_family = scenario$factor_family, 
            rank = q, 
            weights = scenario$factor_weights,
            algo = algo,
            start = init,
            em_control = control,
            ngq = 10,
            eval_likeli = TRUE)
      } else {
        adam_control <- adam.control(max_epoch = 3, batch_size = 16, step_size = 0.05)
        sample_control <- sample.control(sample_size = 50, eval_size = 100)
        
        efm(truth$X/scenario$factor_weights, 
            factor_family = scenario$factor_family, 
            rank = q, 
            weights = scenario$factor_weights,
            algo = algo,
            start = init,
            adam_control = adam_control,
            sample_control = sample_control,
            eval_likeli = TRUE)
      }
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n\n")
      return(NULL)
    })
    
    end_time <- Sys.time()
    elapsed <- end_time - start_time
    
    if (!is.null(res)) {
      cat("  RESULTS:\n")
      cat("    Time:", round(as.numeric(elapsed), 2), "seconds\n")
      cat("    Final loss:", tail(res$like_list, 1), "\n")
      cat("    Iterations:", length(res$like_list), "\n")
      cat("  ✅ SUCCESS\n")
    } else {
      cat("  ❌ FAILED\n")
    }
    cat("\n")
  }
  cat("\n")
}

cat("=== DEMONSTRATION COMPLETE ===\n")
