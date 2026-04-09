#!/bin/bash
                          # Memory (adjust as needed, e.g., 4G for 4 GB)
#SBATCH --time=00:10:00                   # Max run time (hh:mm:ss, adjust as needed)
#SBATCH --account=jvmorr1

# Load any necessary modules here

# Run your bash script
#chrom="$1"                # e.g., from Snakemake wildcards
#gwas_info_file="$2"       # path to info file
#ref_bim="$3"              # reference file to use for snp harmonization
#af_thresh="$4"            # allele freq threshold
#sample_size_tol="$5"      # sample size tolerance
#out="$6"                  # output file

#bash bash/1_combine_and_format.sh {wildcards.chrom} {input.gwas_info} {input.ref_bim} {params.af_thresh} {params.sample_size_tol} {output.out}

# need to put into snakemake call.  or maybe not if R is in their conda env

#    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/C100001554_And_Friends_3Metabolites.csv \
bash 1_gather_snps.sh \
    1 \
    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/5e5Sig_Herit_Mets_8ForLDSCStrip.csv \
    /nfs/turbo/sph-jvmorr/GFA_metabolites_2025/gfa_pipeline/bash/keep8_0_get_ss_bounds_20260330_092135/ss_table.tsv  \
    0.05 \
    0.1 \
    final_pass_snps.txt


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
