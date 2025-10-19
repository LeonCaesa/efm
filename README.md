# Computational Approaches for Exponential-Family Factor Analysis {https://arxiv.org/abs/2403.14925}

## R package installation guide

### Core Package Dependencies
The EFM package has **minimal dependencies** and does **NOT require tidyverse**. The core package only depends on:
- MASS (for matrix operations)
- matrixStats (for numerical stability)  
- stats, methods (base R)

Once you have `devtools` installed, you can install the required packages using the following commands:

```r
  devtools::install_github("LeonCaesa/efm", dependencies = TRUE)
  devtools::install_github("carvalho-research/gaussquadr", dependencies = TRUE)
```

### Quick Start (No Tidyverse Required)
```r
# Install core dependencies
Rscript install_dependencies.R

# Run main example (demonstrates all core functionality)
Rscript example.R
```

**Important**: The main EFM functionality and `R/example.R` work entirely with base R and minimal dependencies. Tidyverse packages are only required for reproducing specific case studies (Sections 4.3-4.4) that involve complex data processing and visualization of computer vision and network datasets.


## Reproducing our experiments

The main code for the Exponential Factor Model (EFM) is in `./R/efm.R`, which covers the model and its corresponding optimization algorithms (Sections 1 - 3) of the paper.

### Core Functionality
- `./example.R` - Comprehensive demonstration of all EFM algorithms across exponential families

### Paper Section Reproductions
To reproduce the experimental results (Section 4. Examples and Results) of the paper:
- `./demo/simustudy/` - Section 4.1: Simulated data and optimization efficiency
- `./demo/covstudy/` - Section 4.2: Covariance modeling and simulation  
- `./demo/casestudy/` - Sections 4.3-4.4: Computer vision and multiplex network experiments

Each demo subfolder contains its own README with specific instructions.

## Repository Structure

```
efm/
├── example.R                    # Main demo (NO tidyverse, all algorithms)
├── install_dependencies.R      # One-command setup for reviewers
├── README.md                   # Emphasizes minimal dependencies
├── R/
│   ├── efm.R                   # Core package (enhanced docs, reduced deps)
│   └── utils.R                 # Utilities (custom rmvnorm, fixed bugs)
└── demo/                       # Paper experiment reproductions
    ├── README.md               # Concise high-level overview
    ├── simustudy/              # Section 4.1 (optimization efficiency)
    ├── covstudy/               # Section 4.2 (covariance modeling)
    │   └── submit_cov.sh       # Batch execution script
    └── casestudy/              # Sections 4.3-4.4 (vision & networks)
        ├── README.md           # Clear usage instructions
        └── data/               # ORL face datasets
```
