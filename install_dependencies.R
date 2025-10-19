# Installation script for EFM package dependencies
# Run this script first before running any demo files

cat("Installing required packages for EFM package...\n")

# List of required packages
required_packages <- c(
  "devtools",
  "MASS",
  "matrixStats",
  "methods",
  "stats"
)

# Optional packages (will warn if not available)
optional_packages <- c(
  "testthat",
  "knitr", 
  "rmarkdown"
)

# GitHub packages
github_packages <- list(
  "gaussquadr" = "carvalho-research/gaussquadr"
)

# Function to install packages if not already installed
install_if_missing <- function(packages, from_cran = TRUE) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      if (from_cran) {
        install.packages(pkg, repos = "https://cran.r-project.org")
      }
    } else {
      cat(pkg, "already installed.\n")
    }
  }
}

# Install CRAN packages
cat("\n=== Installing CRAN packages ===\n")
install_if_missing(required_packages)

cat("\n=== Installing optional packages ===\n")
install_if_missing(optional_packages)

# Install GitHub packages
cat("\n=== Installing GitHub packages ===\n")
if (requireNamespace("devtools", quietly = TRUE)) {
  for (pkg_name in names(github_packages)) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      cat("Installing", pkg_name, "from GitHub...\n")
      devtools::install_github(github_packages[[pkg_name]], quiet = TRUE)
    } else {
      cat(pkg_name, "already installed.\n")
    }
  }
} else {
  cat("devtools not available - skipping GitHub packages\n")
}

cat("\n=== Installation Summary ===\n")
cat("Required packages:\n")
for (pkg in required_packages) {
  status <- if (requireNamespace(pkg, quietly = TRUE)) "✓ OK" else "✗ MISSING"
  cat(" ", pkg, ":", status, "\n")
}

cat("\nGitHub packages:\n")
for (pkg_name in names(github_packages)) {
  status <- if (requireNamespace(pkg_name, quietly = TRUE)) "✓ OK" else "✗ MISSING"
  cat(" ", pkg_name, ":", status, "\n")
}

cat("\nInstallation complete! You can now run the demo files.\n")
