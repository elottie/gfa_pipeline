#!/bin/bash
                          # Memory (adjust as needed, e.g., 4G for 4 GB)
#SBATCH --time=00:10:00                   # Max run time (hh:mm:ss, adjust as needed)
#SBATCH --account=jvmorr1

mkdir -p ../gfa_data

#    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/C100001554_And_Friends_3Metabolites.csv \
bash 1_gather_snps.sh \
    1 \
    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/5e5Sig_Herit_Mets_8ForLDSCStrip.csv \
    ../gfa_data/5e5Sig_Herit_Mets8_sample_size_table.tsv \
    0.05 \
    ../gfa_data/5e5Sig_Herit_Mets8_snps_chr1.tsv


#c                  <- as.numeric(args[1])  # chromosome
#gwas_info_file     <- args[2]
#af_thresh          <- as.numeric(args[3])
#sample_size_tol    <- as.numeric(args[4])
#out                <- args[5]

#conda activate snakemake9
#Rscript 1_mod_combine_and_format.R \
#    1 \
#    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/5e5Sig_Herit_Mets_8ForLDSCStrip.csv \
#    0.05 \
#    0.1 \
#    8mets_1_comb_form_output_R.tsv
