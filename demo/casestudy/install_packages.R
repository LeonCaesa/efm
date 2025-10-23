# Installation script for Computer Vision Analysis Case Study
# Run this script before running computer_vision_analysis.R

cat("Installing packages required for Computer Vision Analysis case study...\n")

# Core packages needed for the analysis
required_packages <- c(
  "devtools",
  "glmnet",
  "MASS",
  "R.matlab",
  "tidyverse",
  "mvtnorm"
)

# Visualization and analysis packages
analysis_packages <- c(
  "plotly",
  "Rtsne",
  "fpc",
  "nnet",
  "rpart", 
  "caret",
  "HandTill2001",
  "gridExtra"
)

# Network analysis packages (for multiplex_network_analysis.R)
network_packages <- c(
  "multinet",
  "dplyr",
  "ggplot2", 
  "stringr",
  "fields",
  "latex2exp"
)

# Special handling for igraph compatibility
special_packages <- list(
  "igraph" = "1.5.1"  # Use older compatible version
)

# Optional packages (will warn if not available but won't stop execution)
optional_packages <- c(
  "snedata"
)

# GitHub packages
github_packages <- list(
  "dmf" = "carvalho-research/dmf"
)

# Function to install packages if not already installed
install_if_missing <- function(packages, from_cran = TRUE, package_type = "required") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      if (from_cran) {
        tryCatch({
          install.packages(pkg, repos = "https://cran.r-project.org")
          cat("✓", pkg, "installed successfully\n")
        }, error = function(e) {
          if (package_type == "optional") {
            cat("⚠ Warning: Could not install optional package", pkg, "\n")
          } else {
            cat("✗ Error installing", pkg, ":", e$message, "\n")
          }
        })
      }
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Install core packages
cat("\n=== Installing core packages ===\n")
install_if_missing(required_packages, package_type = "required")

cat("\n=== Installing analysis packages ===\n")
install_if_missing(analysis_packages, package_type = "required")

cat("\n=== Installing network analysis packages ===\n")
install_if_missing(network_packages, package_type = "required")

cat("\n=== Installing compatibility packages ===\n")
for (pkg_name in names(special_packages)) {
  version <- special_packages[[pkg_name]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat("Installing", pkg_name, "version", version, "for compatibility...\n")
    tryCatch({
      # Try to install specific version
      devtools::install_version(pkg_name, version = version, repos = "https://cran.r-project.org")
      cat("✓", pkg_name, "version", version, "installed successfully\n")
    }, error = function(e) {
      cat("⚠ Could not install specific version, trying latest...\n")
      install.packages(pkg_name, repos = "https://cran.r-project.org")
    })
  } else {
    cat("✓", pkg_name, "already installed\n")
  }
}

cat("\n=== Installing optional packages ===\n")
install_if_missing(optional_packages, package_type = "optional")

# Install GitHub packages
cat("\n=== Installing GitHub packages ===\n")
if (requireNamespace("devtools", quietly = TRUE)) {
  for (pkg_name in names(github_packages)) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      cat("Installing", pkg_name, "from GitHub...\n")
      tryCatch({
        devtools::install_github(github_packages[[pkg_name]], quiet = TRUE)
        cat("✓", pkg_name, "installed successfully from GitHub\n")
      }, error = function(e) {
        cat("✗ Error installing", pkg_name, "from GitHub:", e$message, "\n")
        cat("You can try installing manually with: devtools::install_github('", github_packages[[pkg_name]], "')\n", sep = "")
      })
    } else {
      cat("✓", pkg_name, "already installed\n")
    }
  }
} else {
  cat("✗ devtools not available - cannot install GitHub packages\n")
  cat("Please install devtools first: install.packages('devtools')\n")
}

# Installation summary
cat("\n=== Installation Summary ===\n")
all_packages <- c(required_packages, analysis_packages, network_packages, names(special_packages), optional_packages, names(github_packages))

installed_count <- 0
missing_count <- 0

for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
    installed_count <- installed_count + 1
  } else {
    cat("✗", pkg, "(missing)\n")
    missing_count <- missing_count + 1
  }
}

cat("\nSummary:", installed_count, "packages installed,", missing_count, "packages missing\n")

if (missing_count == 0) {
  cat("\n🎉 All packages installed successfully!\n")
  cat("You can now run computer_vision_analysis.R\n")
} else {
  cat("\n⚠ Some packages are missing. The analysis may not work properly.\n")
  cat("Please check the error messages above and install missing packages manually.\n")
}

# Special note about dmf package
if (!requireNamespace("dmf", quietly = TRUE)) {
  cat("\n📝 Note: The 'dmf' package is critical for this analysis.\n")
  cat("If automatic installation failed, try:\n")
  cat("  devtools::install_github('carvalho-research/dmf')\n")
}
