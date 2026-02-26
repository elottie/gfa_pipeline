#!/bin/bash
#SBATCH --mem=8G                          # Memory (adjust as needed, e.g., 4G for 4 GB)
#SBATCH --time=00:20:00                   # Max run time (hh:mm:ss, adjust as needed)
#SBATCH --account=jvmorr1

# Load any necessary modules here

# Run your bash script
# bash format_flat_chrom.sh file \
#   chrom af_thresh \
#   snp_name pos_name chrom_name A1_name A2_name beta_hat_name se_name p_value_name af_name sample_size_name \
#   effect_is_or \
#   output.txt

bash format_flat_chrom.sh /nfs/turbo/sph-jvmorr/gwas_summary_statistics/METSIM/with_rsid/C100001554/C100001554_regenie_rsid.tsv.gz \
    1 0.05 \
    RSID GENPOS CHROM ALLELE1 ALLELE0 BETA SE LOG10P A1FREQ N \
    FALSE \
    C100001554_format_flat_chrom.tsv

# in R

#conda activate snakemake9
#Rscript run_format_flat_chrom.R
