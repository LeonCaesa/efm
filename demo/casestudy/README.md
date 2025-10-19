# EFM Case Studies (Paper Sections 4.3-4.4)

This directory contains case study demonstrations for computer vision and network analysis applications.

## Essential Files

### Section 4.3: Computer Vision
- **`computer_vision_analysis.R`** - Computer vision experiments using ORL face dataset
- **`data/`** - Contains ORL face datasets (32x32 and 64x64 pixels)

### Section 4.4: Multiplex Network Analysis  
- **`multiplex_network_analysis.R`** - Multiplex network analysis using AUCS dataset
- **`utilities.R`** - Utility functions used by the main scripts

## Usage

```bash
cd casestudy/
Rscript computer_vision_analysis.R    # Computer vision experiment
Rscript multiplex_network_analysis.R  # Network analysis experiment
```

## Notes

All essential files for reproducing paper sections 4.3-4.4 are included. Additional visualization scripts have been removed to maintain a clean, publication-ready repository.

## Output

Results and plots are saved to local `results/` directory within this folder.

## Dependencies

These scripts require additional packages beyond the core EFM dependencies:
- Computer vision: R.matlab, tidyverse, plotly, etc.
- Network analysis: multinet, igraph, plotly, etc.

All dependencies are automatically installed when the scripts run.
