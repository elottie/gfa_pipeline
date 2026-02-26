#!/bin/bash

echo "in format_flat_chrom.sh:"

set -x

source get_col.sh
source gwas_format.sh

# non great lakes people will need to modify
#module load Bioinformatics
#module load bcftools

# Usage:
# bash format_flat_chrom.sh file chrom af_thresh snp_name pos_name chrom_name \
#   A1_name A2_name beta_hat_name se_name p_value_name af_name sample_size_name effect_is_or output.txt
# bash format_flat_chrom.sh data.tsv.gz 1 0.05 SNP POS CHR A1 A2 BETA SE PVALUE AF SS FALSE output.txt

echo " - loaded modules"

set -euo pipefail
trap 'echo "script failed with exit code $? at line $LINENO" >&2' ERR

echo " - set pipefail"

file="$1"
chrom="$2"
af_thresh="$3"
snp_name="$4"
pos_name="$5"
chrom_name="$6"
A1_name="$7"
A2_name="$8"
beta_hat_name="$9"
se_name="${10}"
p_value_name="${11}"
af_name="${12}"
sample_size_name="${13}"
effect_is_or="${14}"   # "TRUE" or "FALSE"
output="${15}"

echo " - read function inputs"
parse_header "$file"

snp_col_ind=$(get_col "$snp_name")
pos_col_ind=$(get_col "$pos_name")
chrom_col_ind=$(get_col "$chrom_name")
A1_col_ind=$(get_col "$A1_name")
A2_col_ind=$(get_col "$A2_name")
beta_col_ind=$(get_col "$beta_hat_name")
se_col_ind=$(get_col "$se_name")
p_col_ind=$(get_col "$p_value_name")
af_col_ind=$(get_col "$af_name")
ss_col_ind=$(get_col "$sample_size_name")

echo " - parsed gwas file header"

{
  echo "Mem used after step 1: parsing gwas file header"
  awk '/VmRSS/{print}' /proc/self/status
} >> format_flat_chrom_sh_memory.log

# step 2:  read in whole files, filter to chrom, and deduplicate
# best to do upfront to reuse
# written here instead of above w/ headers for extra clarity at the cost of a little efficiency

# read in whole gwas file, filter to chrom and deduplicate (rm all rsids which appear >1x)
# this code keeps header in
# HARDCODING OF DELIMITER
#delim=$(gzip -cd "$file" | head -n 1 | grep -o $'\t' | wc -l)
#if (( delim > 0 )); then
#  awk_delim='\t'
#elif gzip -cd "$file" | head -n 1 | grep -q ','; then
#  awk_delim=','
#else
#  awk_delim=' '
#fi

#Then use:

#awk -F"$awk_delim" ...

# have to make temp files for awk
file_filt=$(mktemp)

# unzip, filter, and deduplicate gwas file
if [[ "$file" == *.gz ]]; then
    gzip -cd "$file" \
    | awk -F'\t' -v col="$chrom_col_ind" -v chrom="$chrom" 'NR==1 || $col == chrom' \
    | awk -F'\t' -v snp_col="$snp_col_ind" '
            NR==1 {print; next}
            {count[$snp_col]++; lines[NR]=$0; keys[NR]=$snp_col}
            END {for (i=2; i<=NR; i++) if (count[keys[i]]==1) print lines[i]}
      ' > "$file_filt"
else
    awk -F'\t' -v col="$chrom_col_ind" -v chrom="$chrom" 'NR==1 || $col == chrom' "$file" \
    | awk -F'\t' -v snp_col="$snp_col_ind" '
       NR==1 {print; next}
       {count[$snp_col]++; lines[NR]=$0; keys[NR]=$snp_col}
       END {for (i=2; i<=NR; i++) if (count[keys[i]]==1) print lines[i]}
      ' > "$file_filt"
fi

{
  echo "Mem used after step 2b: read and filter gwas file"
  awk '/VmRSS/{print}' /proc/self/status
} >> format_flat_chrom_sh_memory.log

# NEED MORE MISSING HANDLING

# A1 is effect allele!!!

# step 3: harmonization of gwas file to reference file

# first, create output file:
# Get parent directory name
parentdir=$(basename "$(dirname "$file")")        # C100001554
# Get base file name without extensions
fname=$(basename "$file")                         # C100001554_regenie_rsid.tsv.gz
# removes everything after first dot, may be too agressive!!!
fname_noext="${fname%%.*}"                        # C100001554_regenie_rsid
# Compose new output file name
harmon="${fname_noext}_chrom${chrom}_harmon_snps.tsv"


# Filter rows by allele frequency
awk -F'\t' -v af_col="$af_col_ind" -v af_thresh="$af_thresh" \
    'NR==1 {print; next} ($af_col > af_thresh && $af_col < (1 - af_thresh)) {print}' \
    "$file_filt" > "${file_filt}_af_filt"

# If effect is odds ratio, add log transformed column "beta"
if [ "${effect_is_or,,}" = "true" ]; then
    awk -F'\t' -v beta_col="$beta_col_ind" 'NR==1 {print $0 "\tbeta"; next} {print $0 "\t" log($beta_col)}' \
        "${file_filt}_af_filt" > "${file_filt}_final"
    beta_hat="beta"
else
    cp "${file_filt}_af_filt" "${file_filt}_final"
    beta_hat="$beta_hat_name"
fi

head "${file_filt}_final"

# now call gwas_format
# Usage: gwas_format <input> <output> <snp> <beta_hat> <se> <A1> <A2> [chrom] [pos] [p_value] [sample_size] [allele_freq] [compute_pval]

#snp_col_ind=$(get_col "$snp_name")
#pos_col_ind=$(get_col "$pos_name")
#chrom_col_ind=$(get_col "$chrom_name")
#A1_col_ind=$(get_col "$A1_name")
#A2_col_ind=$(get_col "$A2_name")
#beta_col_ind=$(get_col "$beta_hat_name")
#se_col_ind=$(get_col "$se_name")
#p_col_ind=$(get_col "$p_value_name")
#af_col_ind=$(get_col "$af_name")
#ss_col_ind=$(get_col "$sample_size_name")

gwas_format "${file_filt}_final" "$output" "$snp_name" "$beta_hat_name" "$se_name" "$A1_name" "$A2_name" "$chrom_name" "$pos_name" "$p_value_name" \
"$sample_size_name" "$af_name" "TRUE"

echo " - formatted output written to $output"
