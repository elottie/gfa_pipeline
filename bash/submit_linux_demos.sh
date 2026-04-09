#!/bin/bash
#SBATCH --job-name=filter_chr
#SBATCH --output=filter_chr.%j.out
#SBATCH --error=filter_chr.%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --account=jvmorr1

bash linux_demos.sh
