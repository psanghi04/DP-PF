# DP-PF

## Project Overview  
This project implements Bayesian linear regression and MCMC workflows, real-data analyses with privacy-preserving methods, and Approximate Bayesian Computation via rejection sampling. Four modules are provided: Linear Regression, MCMC, Real Data Analysis, and Rejection ABC.

## Dependencies  
- R (>= 4.0)  
- Required R packages: `...`  
- SLURM scheduler for batch scripts

## Linear Regression  
Scripts in `linear_regression/`.

### lr_data_org.R  
This script organizes and post-processes data from linear regression experiments. It ensures the data is formatted correctly for downstream analysis.  
**Usage:**  
Rscript lr_data_org.R

### lr_slurm.sh  
A SLURM batch script for submitting linear regression experiments to a high-performance computing cluster. It automates the execution of R scripts for both conjugate and non-conjugate prior models.  
**Usage:**  
sbatch lr_slurm.sh

### conjugate_prior/lr_conj.R  
<Description> 

**Usage:**  
Rscript conjugate_prior/lr_conj.R

### non_conjugate_prior/lr_nonconj.R  
<Description>  

**Usage:**  
Rscript non_conjugate_prior/lr_nonconj.R

## MCMC  
Scripts in `mcmc/`.

### mcmc_data_org.R  
<Description>  

**Usage:**  
Rscript mcmc_data_org.R

### mcmc_slurm.sh  
A SLURM batch script for submitting MCMC experiments to a high-performance computing cluster. It automates the execution of R scripts for both conjugate and non-conjugate prior models.  
**Usage:**  
sbatch mcmc_slurm.sh

### conjugate_prior/mcmc_conj.R  
<Description>

**Usage:**  
Rscript conjugate_prior/mcmc_conj.R  
[Note: uses `...` for sampling]

### non_conjugate_prior/mcmc_nonconj.R  
<Description>  

**Usage:**  
Rscript non_conjugate_prior/mcmc_nonconj.R

## Real Data Analysis  
Scripts in `real_data_analysis/`.

### Canadian_Census_DP_PF_ObjectivePerturb_Knorm.R  
Analyzes real-world Canadian census data using the DP-PF framework with objective perturbation and K-norm regularization. This script applies privacy-preserving Bayesian inference techniques.  
**Usage:**  
Rscript Canadian_Census_DP_PF_ObjectivePerturb_Knorm.R  
[Input: reads census CSV from `data/real/canadian_census.csv`]  
[Output: writes analysis results to `...`]

### plotter.R  
Generates visualizations for the results of real data analysis. This script produces plots such as posterior distributions and parameter estimates.  
**Usage:**  
Rscript plotter.R

### real_data_slurm.sh  
A SLURM batch script for submitting real data analysis tasks to a high-performance computing cluster. 
**Usage:**  
sbatch real_data_slurm.sh

## Rejection ABC  
Scripts in `rejection_abc/`.

### rejection_abc.R  
<Description>  

**Usage:**  
Rscript rejection_abc.R

### rejection_data_org.R  
<Description>  

**Usage:**  
Rscript rejection_data_org.R

### rejection_slurm.sh  
A SLURM batch script for submitting rejection ABC experiments to a high-performance computing cluster.  
**Usage:**  
sbatch rejection_slurm.sh

## How to Cite  
If you use this code in published work, please cite:  
> ...

## Contact  
For questions or issues, please open an issue on the repository or contact ...
