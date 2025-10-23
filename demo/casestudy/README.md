# EFM Case Studies (Paper Sections 4.3-4.4)

This directory contains case study demonstrations for computer vision and network analysis applications.

## Essential Files

### Dependencies
- **`install_packages.R`** - Installation script for all required packages (run this first!)

### Section 4.3: Computer Vision
- **`computer_vision_analysis.R`** - Computer vision experiments using ORL face dataset
- **`data/`** - Contains ORL face datasets (32x32 and 64x64 pixels)

### Section 4.4: Multiplex Network Analysis  
- **`multiplex_network_analysis.R`** - Multiplex network analysis using AUCS dataset
- **`multiplex_network_analysis_simple.R`** - Working alternative for network analysis
- **`utilities.R`** - Utility functions used by the main scripts

## Quick Start

**IMPORTANT**: Before running the analysis scripts, install all required dependencies:

```bash
cd casestudy/
Rscript install_packages.R           # Install all required packages first
```

Then run the analysis:

```bash
Rscript computer_vision_analysis.R           # Computer vision experiment
Rscript multiplex_network_analysis_simple.R  # Network analysis experiment (simplified version)
```

**Note**: The original `multiplex_network_analysis.R` has compatibility issues with newer versions of the `igraph` package. We provide `multiplex_network_analysis_simple.R` as a working alternative that demonstrates the same EFM methodology using synthetic data.

## Dependencies

### Critical Dependencies
The computer vision analysis requires the `dmf` package from GitHub, which is automatically installed by the `install_packages.R` script. This package is essential for the analysis to work properly.

### Complete Package List
- **Core packages**: devtools, glmnet, MASS, R.matlab, tidyverse, mvtnorm
- **Analysis packages**: plotly, Rtsne, fpc, nnet, rpart, caret, HandTill2001, gridExtra
- **GitHub packages**: dmf (from carvalho-research/dmf)
- **Optional**: snedata (analysis will work without this)

### Manual Installation (if needed)
If the automatic installation fails, you can install the critical `dmf` package manually:

```r
devtools::install_github('carvalho-research/dmf')
```

## Notes

All essential files for reproducing paper sections 4.3-4.4 are included. The `install_packages.R` script ensures all dependencies are properly installed, including the critical `dmf` package that was previously causing errors.

## Output

Results and plots are saved to local files within this directory. The computer vision analysis will generate visualization plots and save model results.
