#!/bin/bash
                          # Memory (adjust as needed, e.g., 4G for 4 GB)
#SBATCH --time=00:10:00                   # Max run time (hh:mm:ss, adjust as needed)
#SBATCH --account=jvmorr1

#gwas_info_file="$1"       # path to info file
#ss_tol="$2"      # sample size tolerance
#out="$3"                  # output file

bash 0_get_ss_bounds.sh \
    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/5e5Sig_Herit_Mets_8ForLDSCStrip_1.csv \
    ss_table.tsv
