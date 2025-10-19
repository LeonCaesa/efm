# Paper Experiment Reproductions

This directory contains scripts to reproduce the experimental results from "Computational Approaches for Exponential-Family Factor Analysis" (arXiv:2403.14925).

## Prerequisites
```bash
Rscript ../install_dependencies.R
```

## Paper Sections

### Section 4.1: Simulated Data (`simustudy/`)
Optimization efficiency comparison across algorithms and families.

### Section 4.2: Covariance Modeling (`covstudy/`) 
Covariance modeling and simulation experiments.

### Section 4.3-4.4: Case Studies (`casestudy/`)
Computer vision (ORL faces) and multiplex network analysis.

## Usage

Each subfolder contains:
- **Scripts** to run experiments
- **README.md** with detailed instructions
- **Local results/** directories for output

**Note**: Core EFM functionality demonstration is now in `../example.R` (no tidyverse required).
