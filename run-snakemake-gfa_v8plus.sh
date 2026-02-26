#!/bin/bash

snakemake --profile snakemake-profile-slurm --executor slurm -s Snakefile_gfa
#snakemake --profile snakemake-profile-slurm -s Snakefile_gfa

#snakemake \
#   -s Snakefile_gfa \
#   --keep-going \
#   --notemp \
#   --jobs 96 \
#   --max-jobs-per-second 5 \
#   --latency-wait 30 \
#   --default-resources mem_mb=5000 runtime=360 account=jvmorr1 \
#   --executor slurm
