#!/bin/bash
#SBATCH --job-name=my_r_job       # Job name
#SBATCH --account=your_account    # Account name (e.g., your PI's shortcode)
#SBATCH --partition=standard      # Partition (sometimes "standard", "gpu", etc.)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --mem=2G                 # Memory (adjust as needed)

#module purge
#module load R
# use conda env

# Optional: Print start time and host
echo "Starting job on $(hostname) at $(date)"

# Run your script
Rscript ldsc_strip.R

echo "Job finished at $(date)"
