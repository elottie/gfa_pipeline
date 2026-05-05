#!/bin/bash
#SBATCH --job-name=make_nice_data       # Job name
#SBATCH --account=jvmorr1    # Account name (e.g., your PI's shortcode)
#SBATCH --partition=standard      # Partition (sometimes "standard", "gpu", etc.)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --mem=3G                 # Memory (adjust as needed)

#module purge
#module load R
# use conda env

# Optional: Print start time and host
echo "Starting job on $(hostname) at $(date)"

# Run your script
Rscript 2_make_nice_data.R

echo "Job finished at $(date)"
