#!/bin/bash
                          # Memory (adjust as needed, e.g., 4G for 4 GB)
#SBATCH --time=00:10:00                   # Max run time (hh:mm:ss, adjust as needed)
#SBATCH --account=jvmorr1

#gwas_info_file="$1"       # path to info file
#ss_tol="$2"      # sample size tolerance
#out="$3"                  # output file

mkdir -p ../gfa_data

bash 0_get_ss_bounds.sh \
    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/First8_Mets_ForLDSCStrip.csv \
    0.1 \
    ../gfa_data/First8_Mets_sample_size_table.tsv
