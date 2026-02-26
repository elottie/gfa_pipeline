#!/bin/bash

# non great lakes people will need to modify
module load Bioinformatics
module load bcftools

# Usage: bash format_ieu_chrom.sh input.vcf.gz chrom af_thresh output.txt
# Usage: bash format_ieu_chrom.sh input.vcf.gz 1 0.05 output.txt

set -euo pipefail

# Input parameters
vcf_file="$1"        # path to .vcf.gz
chrom="$2"           # chromosome value, e.g. "1"
af_thresh="$3"       # allele frequency threshold
output="$4"          # output filename

# Step 1: Extract variants for the chromosome
bcftools view -r "$chrom" "$vcf_file" -Oz -o temp_chr${chrom}.vcf.gz

# Step 2: Convert VCF to tabular format with necessary fields
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%ES\t%SE\t%SS\t%LP\n' temp_chr${chrom}.vcf.gz > temp_chr${chrom}.txt

# Step 3: Filter on allele frequency
awk -v af_thresh="$af_thresh" '{ if ($6 > af_thresh && $6 < (1 - af_thresh)) print }' temp_chr${chrom}.txt > temp_chr${chrom}.filtered.txt

# Step 4: Compute p-value (10^-LP) and format output
awk '{
    lp = $10;
    p_value = 10^(-lp);
    print $3, $7, $8, $5, $4, $1, $2, p_value, $9, $6
}' temp_chr${chrom}.filtered.txt > temp_chr${chrom}.pv.txt

# Step 5: Add header and save (fields: ID ES SE ALT REF CHROM POS p_value SS AF)
echo -e "ID\tES\tSE\tALT\tREF\tCHROM\tPOS\tp_value\tSS\tAF" > "$output"
cat temp_chr${chrom}.pv.txt >> "$output"

# Step 6: Clean up
rm temp_chr${chrom}.vcf.gz temp_chr${chrom}.txt temp_chr${chrom}.filtered.txt temp_chr${chrom}.pv.txt

echo "Output written to $output"
