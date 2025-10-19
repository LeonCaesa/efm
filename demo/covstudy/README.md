# Section 4.2: Covariance Modeling Experiments

This directory contains scripts to reproduce the covariance modeling and simulation experiments from Section 4.2 of the paper, demonstrating EFM's ability to model complex covariance structures across different exponential families.

## Files

- **`covexp.R`** - Main covariance modeling experiment
- **`evalcov.R`** - Covariance evaluation and metrics computation
- **`plotcov.R`** - Visualization of covariance results
- **`replicate_gaussian2008.R`** - Replication of Gaussian 2008 study
- **`submit_cov.sh`** - Batch execution script for cluster computing

## Usage

### Prerequisites
```bash
cd ../../
Rscript install_dependencies.R
```

### Running Experiments

#### Basic Usage
```bash
cd demo/covstudy/
Rscript covexp.R [exp_idx] [n_repeats] [d]
```

#### Parameters
- **`exp_idx`**: Experiment type (exponential family)
  - `"1"` = Negative Binomial (dispersion=20)
  - `"2"` = Quasi-Poisson (variable dispersion)
  - `"3"` = Binomial (with random weights)
  - `"4"` = Poisson (dispersion=1)
- **`n_repeats`**: Repetition number for random seed (1, 2, 3, ...)
- **`d`**: Number of variables (e.g., 16, 50, 100)

#### Example Commands
```bash
# Test Poisson family, repetition 1, d=16
Rscript covexp.R 4 1 16

# Test Quasi-Poisson family, repetition 2, d=50
Rscript covexp.R 2 2 50

# Test Negative Binomial, repetition 1, d=100
Rscript covexp.R 1 1 100
```

### Batch Execution
For multiple runs (cluster computing):
```bash
bash submit_cov.sh
```

### Evaluation and Plotting
```bash
Rscript evalcov.R     # Compute covariance metrics
Rscript plotcov.R     # Create visualizations
```

## Output

Results are saved to `results/[family]/` directories:
- **Truth data**: `truth_[d]_[n_repeats].RData`
- **EM results**: `fagqem_[d]_[n_repeats].RData`

## Experiment Design

### Data Generation Parameters (from Jianqing Fan 2008)
- **Sample size**: n = 756
- **Factors**: q = 3
- **Prior parameters**: Realistic values from literature
- **Families**: Variable dispersion for quasi-families

### Algorithm Configuration
- **EM algorithm**: Gaussian quadrature (15 nodes)
- **Likelihood evaluation**: Enabled for convergence tracking
- **Iterations**: Adaptive based on convergence

## Notes

This experiment demonstrates EFM's ability to model complex covariance structures that arise in high-dimensional data with non-Gaussian distributions, extending traditional factor analysis approaches as discussed in Section 4.2 of the paper.
