#!/bin/bash

# Users must get csvtk:
#conda activate <your-gfa-env>
#conda install csvtk
#csvtk version

set -x

set -euo pipefail

# read input arguments - - -

#bash bash/1_combine_and_format.sh {wildcards.chrom} {input.gwas_info} {input.ref_bim} {params.af_thresh} {params.sample_size_tol} {output.out}

chrom="$1"                # e.g., from Snakemake wildcards
gwas_info_file="$2"       # path to info file
af_thresh="$3"            # allele freq threshold
sample_size_tol="$4"      # sample size tolerance
out="$5"                  # output file

# temp workspace for generated files
workdir="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$workdir"

# just once, need to get rid of lovely windows carriage returns
dos2unix "$gwas_info_file"
# for them to undo it, they can do
#unix2dos "$gwas_info_file"
# most of the time, they don't need to undo it:  Excel, R, Python recognize Unix endings fine.  just an issue for simple things like Notepad

# add the af thing.  if no af in columns, make an na column of it

# read user metadata file - - -

# Get the number of rows/traits
num_traits=$(awk 'NR>1' "$gwas_info_file" | wc -l)

# Get the header line and column indices
source get_col.sh
delimiter=$(get_file_delimiter "$gwas_info_file")
parse_header "$gwas_info_file" "$delimiter"      # build the associative array col_indices
# snp_col=$(get_col "SNP")         # retrieve the index for column named "SNP"
# awk -v snp_col="$snp_col" ...

col_raw_data=$(get_col "raw_data_path")
col_snp=$(get_col "snp")
col_pos=$(get_col "pos")
col_chrom=$(get_col "chrom")
col_A1=$(get_col "A1")
col_A2=$(get_col "A2")
col_beta_hat=$(get_col "beta_hat")
col_se=$(get_col "se")
col_pval=$(get_col "p_value")
col_af=$(get_col "af")
col_sample_size=$(get_col "sample_size")
col_pub_sample_size=$(get_col "pub_sample_size")
col_effect_is_or=$(get_col "effect_is_or")
col_name=$(get_col "name")

#sqlite3 my_gwas.db <<'EOF'
#DROP TABLE IF EXISTS gwas;
#CREATE TABLE gwas (
#  chrom TEXT,
#  snp TEXT,
#  REF TEXT,
#  ALT TEXT,
#  beta REAL,
#  se REAL,
#  sample_size INTEGER,
#  af REAL,
#  z REAL,
#  trait_name TEXT
#);
#EOF

# loop over gwas files and format - - -

trait_files=()
#for ((i=2; i<=num_traits+1; i++)); do
for ((i=2; i<=2; i++)); do
    # Read info fields for trait i (awk column indexing, adjust as needed for TSV)
    f=$(awk -F, -v row="$i" -v col="$col_raw_data" 'NR==row {print $col}' "$gwas_info_file")
    snp=$(awk -F, -v row="$i" -v col="$col_snp" 'NR==row {print $col}' "$gwas_info_file")
#    pos=$(awk -F, -v row="$i" -v col="$col_pos" 'NR==row {print $col}' "$gwas_info_file")
    chrn=$(awk -F, -v row="$i" -v col="$col_chrom" 'NR==row {print $col}' "$gwas_info_file")
    A1=$(awk -F, -v row="$i" -v col="$col_A1" 'NR==row {print $col}' "$gwas_info_file")
    A2=$(awk -F, -v row="$i" -v col="$col_A2" 'NR==row {print $col}' "$gwas_info_file")
#    beta_hat=$(awk -F, -v row="$i" -v col="$col_beta_hat" 'NR==row {print $col}' "$gwas_info_file")
#    se=$(awk -F, -v row="$i" -v col="$col_se" 'NR==row {print $col}' "$gwas_info_file")
#    pval=$(awk -F, -v row="$i" -v col="$col_pval" 'NR==row {print $col}' "$gwas_info_file")
#    af=$(awk -F, -v row="$i" -v col="$col_af" 'NR==row {print $col}' "$gwas_info_file")
#    sample_size=$(awk -F, -v row="$i" -v col="$col_sample_size" 'NR==row {print $col}' "$gwas_info_file")
#    effect_or=$(awk -F, -v row="$i" -v col="$col_effect_is_or" 'NR==row {print tolower($col)}' "$gwas_info_file")
#    pub_sample_size=$(awk -F, -v row="$i" -v col="$col_pub_sample_size" 'NR==row {print $col}' "$gwas_info_file")
    trait_name=$(awk -F, -v row="$i" -v col="$col_name" 'NR==row {print $col}' "$gwas_info_file")

    trait_out="$workdir/${trait_name}.final.tsv"

    delimiter=$(get_file_delimiter "$f")
    echo "got delimiter"
    parse_header "$f" "$delimiter"
    echo "parsed header"
    chr_col=$(get_col "$chrn")
    echo "$chr_col"

    if [[ "$f" == *.vcf.gz || "$f" == *.vcf.bgz ]]; then
        echo "Calling format_ieu_chrom (external): $f $chrom $af_thresh"
        format_ieu_chrom "$f" "$chrom" "$af_thresh" > "$trait_out"
    else
    # STEP 1: Harmonize, chrom filter, fill sample size
    # add elif for .gz
         (
         zcat "$f" \
             | awk -F"$delimiter" -v col="$chr_col" -v chrom_val="$chrom" -v OFS="\t" 'NR==1 || $col == chrom_val' \
             | awk -F"\t" -v OFS="\t" \
#                 -v compute_pval="TRUE" \
                 -v snp_name="$snp" \
#                 -v beta_name="$beta_hat" \
#                 -v se_name="$se" \
                 -v A1_name="$A1" \
                 -v A2_name="$A2" \
#                 -v pos_name="$pos" \
#                 -v pval_name="$pval" \
#                 -v ss_name="$sample_size" \
#                 -v af_name="$af" \
                 -f remove_invalid_variants.awk
        ) > "$workdir/${trait_name}.gwasform.tsv"
        # connect
        bash make_snp_table.sh

             #| awk -F"\t" -v OFS="\t" \
             #    -v chrom="$chrom" \
             #    -v snp_name="$snp" \
             #    -v chrom_name="$chrn" \
             #    -v af_name="$af" \
             #    -v af_thresh="$af_thresh" \
             #    -v effect_is_or="$effect_or" \
             #    -f format_flat_chrom.awk \
             #| awk -v pub_ss="$pub_sample_size" -f fill_sample_size.awk
        #) > "$workdir/${trait_name}.filled.tsv"

        # STEP 2: Determine sample_size column index and median
        #ss_col=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_size") print i}' "$workdir/${trait_name}.filled.tsv")
        #median=$(awk -v ss_col="$ss_col" 'NR>1 && $ss_col!="NA"{a[++N]=$ss_col} END{n=asort(a); m=int(n/2); print a[m]}' "$workdir/${trait_name}.filled.tsv")
        #lower=$(awk -v m="$median" -v tol="$sample_size_tol" 'BEGIN{print (1-tol)*m}')
        #upper=$(awk -v m="$median" -v tol="$sample_size_tol" 'BEGIN{print (1+tol)*m}')

        # STEP 3: Apply sample_size filter and add Z
        #awk -v ss_col="$ss_col" -v lower="$lower" -v upper="$upper" -f sample_size_tol_filter.awk "$workdir/${trait_name}.filled.tsv" \
        #   | awk -f add_zscore.awk \
        #   > "$trait_out"
    fi

    # ---- Add trait name column, ready for SQLite import ----
    #awk -v trait="$trait_name" 'NR==1{print $0 "\ttrait_name"; next} {print $0 "\t" trait}' "$trait_out" > "$workdir/${trait_name}.with_trait.tsv"

    #trait_files+=("$trait_out")

done

# concatenate, long format:
# Use all with_trait.tsv files
#first_file=$(ls $workdir/*.with_trait.tsv | head -n 1)

# Print header from first file
#head -n 1 "$first_file" > all_traits_long.tsv

# Print data from all files (skip first line—header)
#for f in $workdir/*.with_trait.tsv; do
#    tail -n +2 "$f"
#done >> all_traits_long.tsv

#for trait_file in "${trait_files[@]}"; do
#    sqlite3 my_gwas.db <<EOF
#.mode tabs
#.import $trait_file gwas
#EOF
#done

# merge GWAS files NOT IN BASH - - -

# Use csvtk, datamash, or join for multi-file full join if possible. Otherwise, note this as an R step.
#echo "Full-joining all trait files: ${trait_files[@]}"

# HARDCODING DELIMITER
#csvtk join --outer-join --tabs -f chrom,snp,REF,ALT "${trait_files[@]}" > temp.tsv
# Replace empty tabs with NA
# This sed trick works for tab-separated files: replaces consecutive tabs or tabs at line ends
#sed 's/\t\t/\tNA\t/g; s/\t$/\tNA/g' temp.tsv > fulldat.tsv

# remove dup SNPs - - -

# Find and filter out duplicate snps (assuming final merged file 'fulldat.tsv')
#awk '
#NR==1 {for(i=1;i<=NF;i++) if($i=="snp") snp=i; print; next;}
#{
#    if(seen[$snp]++ == 0) print;
#}
#' all_traits_long.tsv > fulldat.nodup.tsv

# remove missing traits - - -

# For each row, check Z columns, count NAs, keep only rows where all Z columns are present:
#awk '
#NR==1 {
#    for(i=1;i<=NF;i++) if($i ~ /\.z$/) zcols[++zcount]=i;
#    print
#}
#NR>1 {
#    miss=0;
#    for(j=1;j<=zcount;j++) if($zcols[j]=="NA") miss++;
#    if(miss == 0) print;
#}
#' fulldat.nodup.tsv > "${out}"


