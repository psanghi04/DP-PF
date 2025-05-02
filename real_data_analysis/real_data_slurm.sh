#!/bin/bash
#SBATCH --job-name=real_data_simulation
#SBATCH --output=real_data_simulation_%A_%a.log
#SBATCH --error=real_data_simulation_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
##SBATCH --array=0-4 # Used for running multiple chains
#
# Load R module, if necessary
module load r

SUB_DIR="real_data_analysis/data"
SCRIPT_DIR="real_data_analysis"

# Move Directory To Submission Folder
cd $SUB_DIR
id=$((SLURM_ARRAY_TASK_ID % 5 + 1))

Rscript "$SCRIPT_DIR/canadian_census_dp_pf_objective_perturb_knorm.R" --shape 6 --scale 4 --sigma1 16 --sigma2 16 --epsilon 0.5 --id $id

