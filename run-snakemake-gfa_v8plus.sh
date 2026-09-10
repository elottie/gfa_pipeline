#!/bin/bash

#snakemake --profile snakemake-profile-slurm -s Snakefile_gfa

# for debugging
snakemake --profile snakemake-profile-slurm -s Snakefile_gfa --allowed-rules all

# for if snakemake run did not complete and gives error next time you try to run
#snakemake --profile snakemake-profile-slurm -s Snakefile_gfa --unlock
