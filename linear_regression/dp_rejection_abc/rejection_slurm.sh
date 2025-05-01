#!/bin/bash
#SBATCH --job-name=rejection_simulation
#SBATCH --output=rejection_simulation_%j.log
#SBATCH --error=rejection_simulation_%j.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# Load R module, if necessary
module load r/4.2

SUB_DIR="linear_regression/rejection_abc/data"

# Move Directory To Submission Folder
cd $SUB_DIR

# Set The Working Directory To The Location Of The R Script
SCRIPT_DIR="linear_regression/rejection_abc"

# Define arrays of epsilon and sample_size values
epsilons=(0.1 0.5 1.0 2.0)
sample_sizes=(100 200 300 400 500)
particles=(100 200 400)

# Loop through each combination of epsilon and sample_size
for epsilon in "${epsilons[@]}"; do
  for sample_size in "${sample_sizes[@]}"; do
    for particle in "${particles[@]}"; do
      # Submit a separate job for each combination of epsilon and sample_size
      sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=rejection_simulation_eps${epsilon}_ss${sample_size}
#SBATCH --output=rejection_simulation_eps${epsilon}_ss${sample_size}_%j.log
#SBATCH --error=rejection_simulation_eps${epsilon}_ss${sample_size}_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100

# Load R module, if necessary
# module load r/4.2

# Move Directory To Submission Folder
cd $SUB_DIR

for iter in {1..100}; do
  Rscript "$SCRIPT_DIR/rejection_lr_abc_laplace.R" -d \$iter -e $epsilon -n $sample_size -p $particle -i \$iter
done
EOT
   done
  done
done

echo "Script Done"
