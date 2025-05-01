#!/bin/bash
#SBATCH --job-name=lr_simulation
#SBATCH --output=lr_simulation_%j.log
#SBATCH --error=lr_simulation_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# Load R module, if necessary
module load r/4.2

SUB_DIR="linear_regression/dp_pf/conjugate_prior/data"

# Move Directory To Submission Folder
cd $SUB_DIR

# Set The Working Directory To The Location Of The R Script
SCRIPT_DIR="linear_regression/dp_pf/conjugate_prior"

# Define arrays of epsilon and sample_size values
epsilons=(0.1 0.5 1.0 2.0)
sample_sizes=(500)
particles=(100 200 300 400 500 1000 1200 1400 1600)

# Loop through each combination of epsilon and sample_size
for epsilon in "${epsilons[@]}"; do
  for sample_size in "${sample_sizes[@]}"; do
    for particle in "${particles[@]}"; do
      for iter in {1..100}; do
      # Submit a separate job for each combination of epsilon and sample_size
      sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=lr_conj_simulation_eps${epsilon}_ss${sample_size}_i${iter}
#SBATCH --output=lr_conj_simulation_eps${epsilon}_ss${sample_size}_i${iter}_%j.log
#SBATCH --error=lr_conj_simulation_eps${epsilon}_ss${sample_size}_i${iter}_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100

# Load R module, if necessary
# module load r/4.2

# Move Directory To Submission Folder
cd $SUB_DIR

Rscript "$SCRIPT_DIR/lr_conj.R" -d $iter -e $epsilon -n $sample_size -p $particle -i $iter

EOT
    done
   done
  done
done

echo "Script Done"
