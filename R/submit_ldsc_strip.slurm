#!/bin/bash
#SBATCH --job-name=ldsc_strip       # Job name
#SBATCH --account=jvmorr1    # Account name (e.g., your PI's shortcode)
#SBATCH --partition=standard      # Partition (sometimes "standard", "gpu", etc.)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks

#module purge
#module load R
# use conda env
#SBATCH --mem=3G                 # Memory (adjust as needed)

# Optional: Print start time and host
echo "Starting job on $(hostname) at $(date)"

# Run your script
Rscript 3_R_ldsc_strip.R
#Rscript test_R_ldsc.R

echo "Job finished at $(date)"
