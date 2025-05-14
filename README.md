# DP-PF

## Project Overview  
This repository implements DP-PF, a novel particle-filter algorithm for private Bayesian inference that leverages sequential Monte Carlo methods under differential privacy constraints. DP-PF produces weighted samples that strongly converge to the true posterior, guaranteeing consistency and a central limit theorem, enabling valid uncertainty quantification via asymptotic confidence intervals.

## Repository Structure

- **`linear_regression/`**  
  Conjugate and non-conjugate linear regression experiments using DP-PF, Data Augmentation MCMC, and Rejection-ABC.  

- **`location_scale/`**  
  Location-scale experiments with DP-PF and Rejection-ABC.  
  
- **`real_data_analysis/`**  
  Application of DP-PF to real-world datasets (e.g., Canadian census).

- **`particle_filter_video_visualization`**
  Demonstrates how particles shift towards posterior distribution.

## Dependencies

- **R ≥ 4.0**  
- **SLURM**: for batch job scheduling on HPC clusters (If availabe)

## Installation

```bash
git clone https://github.com/psanghi04/DP-PF.git
cd DP-PF
Rscript requirements.R
```

## Usage

1. **Linear Regression**  
   Each DP method includes `*_conj.R` script to run conjugate prior setting, `*_nonconj.R` script to run non-conjugate prior setting, and `*_slurm.sh` script to submit the simulations to a HPC cluster.

2. **Location-Scale**  
   Each DP method includes `*_location_scale.R` script to run the setting and `*_slurm.sh` script to submit the simulations to a HPC.

3. **Real-Data Analysis**  
   `canadian_census_dp_pf_objective_perturb_knorm.R` performs DP-PF on the 2021 Canadian Census data.

4. **Post-Processing & Visualization**  
   Each module includes a `*_data_org.R` script to collate RData outputs and a `*_plotter.R` script to generate publication-quality tables and figures.

## Citation

If you use DP-PF, please cite our paper:

```bibtex
@article{chen2025particle,
  title={Particle Filter for Bayesian Inference on Privatized Data},
  author={Chen, Yu-Wei and Sanghi, Pranav and Awan, Jordan},
  journal={arXiv preprint arXiv:2505.00877},
  year={2025}
}
@misc{}
```

## License

This software is released under the **MIT License**. See [LICENSE](LICENSE) for details.
