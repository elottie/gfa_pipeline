#!/bin/bash
#SBATCH --job-name=coll_ldsc_strip       # Job name
#SBATCH --account=jvmorr1    # Account name (e.g., your PI's shortcode)
#SBATCH --partition=standard      # Partition (sometimes "standard", "gpu", etc.)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --mem=2688M                 # Memory (adjust as needed)

#module purge
#module load R
# use conda env

# Optional: Print start time and host
echo "Starting job on $(hostname) at $(date)"

# Run your script
Rscript 3_R_collect_ldsc_strips.R
#Rscript test_R_ldsc.R

echo "Job finished at $(date)"
