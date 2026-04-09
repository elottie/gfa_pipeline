#!/bin/bash

#awk -F"\t" -v OFS="\t" -v chr_idx=2 -v chrom_val=1 '
#        $chr_idx == chrom_val
#        ' C999912117_regenie_rsid.tsv > my_simple_awk_output.tsv

f="C999912117_regenie_rsid.tsv"
delimiter=$'\t'
chrom=1
snp="RSID"
A1="ALLELE1"
A2="ALLELE0"
beta_hat="BETA"
se="SE"
sample_size="N"
af="A1FREQ"

cat "$f" |
    awk -F"$delimiter" -v chr_idx=2 -v chrom_val="$chrom" -v OFS="\t" '
        NR==1 ||  $chr_idx == chrom_val
        ' \
    | awk -F"\t" -v OFS="\t" \
        -v snp_name="$snp" -v A1_name="$A1" -v A2_name="$A2" \
        -v beta_name="$beta_hat" -v se_name="$se" \
        -v ss_name="$sample_size" -v af_name="$af" \
        -f remove_invalid_variants.awk \
    > my_simple_awk_output.tsv
