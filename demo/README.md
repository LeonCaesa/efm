# EFM Package Demonstrations

This directory contains demonstration scripts that reproduce the experimental results from the paper "Computational Approaches for Exponential-Family Factor Analysis" (arXiv:2403.14925).

## Quick Start

1. **Install dependencies first:**
   ```bash
   Rscript ../install_dependencies.R
   ```

2. **Run the main example:**
   ```bash
   Rscript example.R
   ```

## File Organization

### Main Example
- **`example.R`** - Comprehensive demonstration of EFM functionality across all exponential families and optimization algorithms

### Paper Section Reproductions

#### Section 4.1: Simulated Data (`simustudy/`)
- **`optiexp.R`** - Optimization efficiency study comparing algorithms
- **`plotopti.R`** - Plotting optimization results
- **`revision_plot.R`** - Additional visualization for paper revision

**Usage:**
```bash
cd simustudy/
Rscript optiexp.R [family_idx] [algo_idx] [d] [q]
```

Where:
- `family_idx`: 1=Poisson, 2=Binomial, 3=Negative Binomial
- `algo_idx`: 1=Posterior Sampling, 2=SML, 3=Laplacian, 4=EM
- `d`: Number of variables
- `q`: Number of factors

#### Section 4.2: Covariance Modeling (`covstudy/`)
- **`covexp.R`** - Covariance modeling experiments
- **`evalcov.R`** - Covariance evaluation
- **`plotcov.R`** - Covariance plotting
- **`replicate_gaussian2008.R`** - Replication of Gaussian 2008 study

**Usage:**
```bash
cd covstudy/
Rscript covexp.R [exp_idx] [n_repeats] [d]
```

Where:
- `exp_idx`: "1"=Negative Binomial, "2"=Quasi-Poisson, "3"=Binomial, "4"=Poisson
- `n_repeats`: Repetition number for random seed
- `d`: Number of variables

#### Section 4.3-4.4: Case Studies (`casestudy/`)
- **`cv_casestudy.R`** - Computer vision experiments (ORL face dataset)
- **`network_casestudy.R`** - Multiplex network experiments
- **`factor_plot.R`** - Factor visualization
- **`util_casestudy.R`** - Utility functions for case studies
- **`data/`** - Contains ORL face datasets (32x32 and 64x64)

**Usage:**
```bash
cd casestudy/
Rscript cv_casestudy.R
```

### Execution Scripts (`exec/`)
- **`submit_cov.sh`** - Batch submission script for cluster computing

## Output

All scripts save results to local `results/` directories within their respective folders. Results are saved as `.RData` files with descriptive names indicating the algorithm, family, and parameters used.

## Dependencies

All scripts assume you have run `../install_dependencies.R` first to install required packages:
- devtools, MASS, matrixStats, gaussquadr (core dependencies)
- Additional packages for specific case studies (automatically handled)

## Notes for Reviewers

- All hardcoded paths have been removed for portability
- Scripts work in any environment after dependency installation
- Results are saved locally rather than to external cluster directories
- All scripts use the optimized EFM package with reduced dependencies
