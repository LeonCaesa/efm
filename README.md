# Computational Approaches for Exponential-Family Factor Analysis {https://arxiv.org/abs/2403.14925}

The main code for the Exponential Factor Model (EFM) is in `./R/efm.R`, which covers the model and its corresponding optimization algorithms (Sections 1 - 3) of the paper.

To reproduce the experimental results (Section 4. Examples and Results) of the paper:
- `./demo/simustudy/` covers Section 4.1 on simulated data to demonstrate the optimization efficiency of the studied algorithms.
- `./demo/covstudy/` covers Section 4.2 on covariance modeling and simulation.
- `./demo/casestudy/` covers Section 4.3 on computer vision experiments and Section 4.4 on the multiplex network experiment.
