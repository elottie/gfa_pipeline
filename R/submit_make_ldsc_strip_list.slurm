#!/bin/bash
#SBATCH --job-name=make_ldsc_strip       # Job name
#SBATCH --account=jvmorr1    # Account name (e.g., your PI's shortcode)
#SBATCH --partition=standard      # Partition (sometimes "standard", "gpu", etc.)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --mem=1G                 # Memory (adjust as needed)

#module purge
#module load python3.11
# use conda env

# Optional: Print start time and host
echo "Starting job on $(hostname) at $(date)"

# Run your script
#python run_make_ldsc_strip_list.py
Rscript make_ldsc_strip_list_3.R

echo "Job finished at $(date)"
