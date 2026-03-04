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
workdir="workdir_$$"
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
header=$(head -1 "$gwas_info_file")
IFS=$',' read -ra cols <<< "$header"    # Change delimiter if TSV

get_col() {
    local name="$1"
    for i in "${!cols[@]}"; do
        if [[ "${cols[$i]}" == "$name" ]]; then
            echo "$((i+1))"
            return
        fi
    done
    echo "-1"
}

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

# loop over gwas files and format - - -

trait_files=()
for ((i=2; i<=num_traits+1; i++)); do
    # Read info fields for trait i (awk column indexing, adjust as needed for TSV)
    f=$(awk -F, -v row="$i" -v col="$col_raw_data" 'NR==row {print $col}' "$gwas_info_file")
    snp=$(awk -F, -v row="$i" -v col="$col_snp" 'NR==row {print $col}' "$gwas_info_file")
    pos=$(awk -F, -v row="$i" -v col="$col_pos" 'NR==row {print $col}' "$gwas_info_file")
    chrn=$(awk -F, -v row="$i" -v col="$col_chrom" 'NR==row {print $col}' "$gwas_info_file")
    A1=$(awk -F, -v row="$i" -v col="$col_A1" 'NR==row {print $col}' "$gwas_info_file")
    A2=$(awk -F, -v row="$i" -v col="$col_A2" 'NR==row {print $col}' "$gwas_info_file")
    beta_hat=$(awk -F, -v row="$i" -v col="$col_beta_hat" 'NR==row {print $col}' "$gwas_info_file")
    se=$(awk -F, -v row="$i" -v col="$col_se" 'NR==row {print $col}' "$gwas_info_file")
    pval=$(awk -F, -v row="$i" -v col="$col_pval" 'NR==row {print $col}' "$gwas_info_file")
    af=$(awk -F, -v row="$i" -v col="$col_af" 'NR==row {print $col}' "$gwas_info_file")
    sample_size=$(awk -F, -v row="$i" -v col="$col_sample_size" 'NR==row {print $col}' "$gwas_info_file")
    effect_or=$(awk -F, -v row="$i" -v col="$col_effect_is_or" 'NR==row {print tolower($col)}' "$gwas_info_file")
    pub_sample_size=$(awk -F, -v row="$i" -v col="$col_pub_sample_size" 'NR==row {print $col}' "$gwas_info_file")
    trait_name=$(awk -F, -v row="$i" -v col="$col_name" 'NR==row {print $col}' "$gwas_info_file")

    #trait_out="$workdir/trait_${i}.tsv"
    trait_out="$workdir/${trait_name}.formatted.tsv"
    # Detect file type and call appropriate formatter
    if [[ "$f" == *.vcf.gz || "$f" == *.vcf.bgz ]]; then
        echo "Calling format_ieu_chrom (external): $f $chrom $af_thresh"
        # Placeholder for actual format_ieu_chrom call
        format_ieu_chrom "$f" "$chrom" "$af_thresh" > "$trait_out"
    else
       # bash bash/format_flat_chrom.sh "$ref_bim" "$f" "$chrom" "$af_thresh" \
       bash format_flat_chrom.sh "$f" "$chrom" "$af_thresh" \
       "$snp" "$pos" "$chrn" "$A1" "$A2" "$beta_hat" "$se" "$pval" "$af" "$sample_size" "$effect_or" "$trait_out"
    fi

    # vars are not what you think
    col_snp_x="snp"
    col_chrom_x="chrom"
    col_A1_x="A1"
    col_A2_x="A2"
    col_beta_hat_x="beta_hat"
    col_se_x="se"
    col_af_x="allele_freq"
    col_sample_size_x="sample_size"

   # do some filtering
   # Fill sample_size if all NA in trait file (can use awk or Python for better handling)
   # debug:  print header
   head -n 1 "$trait_out"

   echo "$trait_out"
   echo "$workdir/${trait_name}.pub_ss.tsv"
   echo "Lines in trait_out BEFORE awk:"
   wc -l "$trait_out"

   awk -v pub_ss="$pub_sample_size" '
       BEGIN{FS=OFS="\t"; ss_col=-1}
       NR==1 {
           for(i=1;i<=NF;i++)
               if($i=="sample_size") ss_col=i;
           print
       }
       NR>1 {
           if($ss_col == "" || $ss_col == "NA") $ss_col=pub_ss;
           print
       }
   ' "$trait_out" > "$workdir/${trait_name}.pub_ss.tsv"

    # Filter by sample_size_tol if finite
    if [[ "$sample_size_tol" != "NA" && "$sample_size_tol" != "" && "$sample_size_tol" != "0" ]]; then
        # Compute median sample size (Awk method, simplified)
        median=$(awk 'NR>1 && $NF!="NA"{a[NR]=$NF} END{n=asort(a); m=int(n/2); print a[m]}' "$trait_out")
        lower=$(awk -v m="$median" -v tol="$sample_size_tol" 'BEGIN{print (1-tol)*m}')
        upper=$(awk -v m="$median" -v tol="$sample_size_tol" 'BEGIN{print (1+tol)*m}')
        awk -v lower="$lower" -v upper="$upper" '
            NR==1 {print}
            NR>1 {if($NF >= lower && $NF <= upper) print}
        ' "$workdir/${trait_name}.pub_ss.tsv" > "$workdir/${trait_name}.ss_tol.tsv"
    else
       # Just copy pub_ss to ss_tol, no sample size tolerance filtering required
        cp "$workdir/${trait_name}.pub_ss.tsv" "$workdir/${trait_name}.ss_tol.tsv"
    fi


    # just cols we need
    awk -F'\t' -v OFS='\t' \
      -v chrom="$col_chrom_x" \
      -v snp="$col_snp_x" \
      -v A2="$col_A2_x" \
      -v A1="$col_A1_x" \
      -v beta="$col_beta_hat_x" \
      -v se="$col_se_x" \
      -v ss="$col_sample_size_x" \
      -v af="$col_af_x" \
      -v trait="$trait_name" \
    '
    BEGIN { OFS="\t" }
    NR==1 {
       # map column names to indices
        for(i=1; i<=NF; i++) {
            if($i==chrom) c_col=i;
            if($i==snp) s_col=i;
            if($i==A2) ref_col=i;
            if($i==A1) alt_col=i;
            if($i==beta) beta_col=i;
            if($i==se) se_col=i;
            if($i==ss) ss_col=i;
            if($i==af) af_col=i;
        }
       # header: rename columns
        print "chrom", "snp", "REF", "ALT", trait".beta", trait".se", trait".ss", trait".af", trait".z";
        next
    }
    {
        # assign variables
        beta = $beta_col;
        se   = $se_col;
        # check everything up front
        if (beta == "" || se == "" || beta == "NA" || se == "NA" || se == "0" || se == 0) {
            z = "NA";
        } else {
            z = beta / se;
        }
        print $c_col, $s_col, $ref_col, $alt_col, beta, se, $ss_col, $af_col, z;
    }
    ' "$workdir/${trait_name}.ss_tol.tsv" > "$workdir/${trait_name}.cut.tsv"

    trait_files+=("$workdir/${trait_name}.cut.tsv")
done

# merge GWAS files NOT IN BASH - - -

# Use csvtk, datamash, or join for multi-file full join if possible. Otherwise, note this as an R step.
echo "Full-joining all trait files: ${trait_files[@]}"

# HARDCODING DELIMITER
csvtk join --outer-join --tabs -f chrom,snp,REF,ALT "${trait_files[@]}" > temp.tsv
# Replace empty tabs with NA
# This sed trick works for tab-separated files: replaces consecutive tabs or tabs at line ends
sed 's/\t\t/\tNA\t/g; s/\t$/\tNA/g' temp.tsv > fulldat.tsv

# remove dup SNPs - - -

# Find and filter out duplicate snps (assuming final merged file 'fulldat.tsv')
awk '
NR==1 {for(i=1;i<=NF;i++) if($i=="snp") snp=i; print; next;}
{
    if(seen[$snp]++ == 0) print;
}
' fulldat.tsv > fulldat.nodup.tsv

# remove missing traits - - -

# For each row, check Z columns, count NAs, keep only rows where all Z columns are present:
awk '
NR==1 {
    for(i=1;i<=NF;i++) if($i ~ /\.z$/) zcols[++zcount]=i;
    print
}
NR>1 {
    miss=0;
    for(j=1;j<=zcount;j++) if($zcols[j]=="NA") miss++;
    if(miss == 0) print;
}
' fulldat.nodup.tsv > "${out}"


