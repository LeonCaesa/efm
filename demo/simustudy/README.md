# Section 4.1: Simulated Data Experiments

This directory contains scripts to reproduce the optimization efficiency experiments from Section 4.1 of the paper, demonstrating the performance comparison of different EFM algorithms across various exponential families.

## Files

- **`optiexp.R`** - Main optimization efficiency experiment
- **`plotopti.R`** - Plotting and visualization of optimization results
- **`revision_plot.R`** - Additional plots for paper revision

## Usage

### Prerequisites
```bash
cd ../../
Rscript install_dependencies.R
```

### Running Experiments

#### Basic Usage
```bash
cd demo/simustudy/
Rscript optiexp.R [family_idx] [algo_idx] [d] [q]
```

#### Parameters
- **`family_idx`**: Exponential family
  - `1` = Poisson
  - `2` = Binomial  
  - `3` = Negative Binomial
- **`algo_idx`**: Optimization algorithm
  - `1` = Posterior Sampling (PS)
  - `2` = Simulated Maximum Likelihood (SML)
  - `3` = Laplacian Approximation (LAPL)
  - `4` = EM Algorithm
- **`d`**: Number of variables (e.g., 50, 100, 512)
- **`q`**: Number of factors (e.g., 3, 10, 50)

#### Example Commands
```bash
# Test Poisson family with Posterior Sampling, d=50, q=3
Rscript optiexp.R 1 1 50 3

# Test Binomial family with EM algorithm, d=100, q=5  
Rscript optiexp.R 2 4 100 5

# Test Negative Binomial with Laplacian, d=512, q=10
Rscript optiexp.R 3 3 512 10
```

### Plotting Results
After running experiments, create plots:
```bash
Rscript plotopti.R
```

## Output

Results are saved to `results/` directory as `.RData` files with naming convention:
```
[algorithm]_[family]_s[sample_size]_d[dimensions]_q[factors]_T[epochs].RData
```

Example: `ps_poisson_s50_d50_q3_T25.RData`

## Experiment Design

- **Sample sizes**: 50, 300, 500 (for PS and SML algorithms)
- **Epochs**: 25 (configurable)
- **Batch size**: 128
- **Evaluation**: Likelihood tracking enabled for convergence analysis

## Notes

This experiment demonstrates the computational efficiency and convergence properties of different optimization algorithms across various exponential family distributions, as discussed in Section 4.1 of the paper.
