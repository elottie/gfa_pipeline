#!/bin/bash

#set -x

set -euo pipefail

# make all awks, joins, and sorts consistent
export LC_ALL=C

# temp workspace for generated files
workdir="0_get_ss_bounds_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$workdir"

memsnap() {
  # prints MaxRSS and AveRSS for the current job step
  sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxRSS,AveRSS -n 2>/dev/null \
    | awk '{print strftime("[%F %T]"), $0}'
}

( while :; do
    date
    ps -u "$USER" -o pid,ppid,rss,cmd --sort=-rss | head -20
    echo "----"
    sleep 2
  done ) > "$workdir/0_ps_rss.log" &
sampler=$!

# read input arguments - - -

#bash bash/1_combine_and_format.sh {wildcards.chrom} {input.gwas_info} {input.ref_bim} {params.af_thresh} {params.sample_size_tol} {output.out}

gwas_info_file="$1"       # path to info file
ss_tol="$2"      # sample size tolerance
out="$3"                  # output file

# just once, need to get rid of lovely windows carriage returns
dos2unix "$gwas_info_file"
# for them to undo it, they can do
#unix2dos "$gwas_info_file"
# most of the time, they don't need to undo it:  Excel, R, Python recognize Unix endings fine.  just an issue for simple things like Notepad

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
col_chrom=$(get_col "chrom")
col_A1=$(get_col "A1")
col_A2=$(get_col "A2")
col_beta_hat=$(get_col "beta_hat")
col_se=$(get_col "se")
col_af=$(get_col "af")
col_sample_size=$(get_col "sample_size")
#col_pub_sample_size=$(get_col "pub_sample_size")
#col_effect_is_or=$(get_col "effect_is_or")
col_name=$(get_col "name")

echo "after extracting col names"
memsnap

# loop over gwas files and format - - -

trait_summary_table="$workdir/$out"
 
# Write header to summary table (only once)
echo -e "trait\tss_low\tss_median\tss_high" > "$trait_summary_table" 

for ((i=2; i<=num_traits+1; i++)); do
#for ((i=2; i<=2; i++)); do
#for ((i=2; i<=3; i++)); do
    # Read info fields for trait i (awk column indexing, adjust as needed for TSV)
    f=$(awk -F, -v row="$i" -v col="$col_raw_data" 'NR==row {print $col}' "$gwas_info_file")
    snp=$(awk -F, -v row="$i" -v col="$col_snp" 'NR==row {print $col}' "$gwas_info_file")
    chrn=$(awk -F, -v row="$i" -v col="$col_chrom" 'NR==row {print $col}' "$gwas_info_file")
    A1=$(awk -F, -v row="$i" -v col="$col_A1" 'NR==row {print $col}' "$gwas_info_file")
    A2=$(awk -F, -v row="$i" -v col="$col_A2" 'NR==row {print $col}' "$gwas_info_file")
    beta_hat=$(awk -F, -v row="$i" -v col="$col_beta_hat" 'NR==row {print $col}' "$gwas_info_file")
    se=$(awk -F, -v row="$i" -v col="$col_se" 'NR==row {print $col}' "$gwas_info_file")
    af=$(awk -F, -v row="$i" -v col="$col_af" 'NR==row {print $col}' "$gwas_info_file")
    sample_size=$(awk -F, -v row="$i" -v col="$col_sample_size" 'NR==row {print $col}' "$gwas_info_file")
#    effect_or=$(awk -F, -v row="$i" -v col="$col_effect_is_or" 'NR==row {print tolower($col)}' "$gwas_info_file")
#    pub_sample_size=$(awk -F, -v row="$i" -v col="$col_pub_sample_size" 'NR==row {print $col}' "$gwas_info_file")
    trait_name=$(awk -F, -v row="$i" -v col="$col_name" 'NR==row {print $col}' "$gwas_info_file")
    
    delimiter=$(get_file_delimiter "$f")

    #trait_out="$workdir/${trait_name}.final.tsv"

    if [[ "$f" == *.vcf.gz || "$f" == *.vcf.bgz ]]; then
        echo "Calling format_ieu_chrom (external): $f $chrom $af_thresh"
    #    format_ieu_chrom "$f" "$chrom" "$af_thresh" > "$trait_out"

    else
        echo "maxrss before using make_filt_data"
        memsnap

        make_ss_data() {
            zcat "$f" |
                awk -F"$delimiter" -v OFS="\t" \
                        -v snp_name="$snp" -v A1_name="$A1" -v A2_name="$A2" \
                        -v beta_name="$beta_hat" -v se_name="$se" \
                        -v ss_name="$sample_size" -v af_name="$af" \
                        -f remove_invalid_variants.awk \
		| awk -F"\t" -v OFS="\t" 'NR>1 {print $1,$6}' \
                | sort -T "$workdir" -S 200M -t $'\t' -k1,1 -u \
		| awk -F"\t" -v OFS="\t" '$2 ~ /^[0-9.]+$/ {print $2}' \
		| sort -T "$workdir" -S 200M -n
	}

        # 1. Sample size statistics for trait summary table

        #make_filt_data | awk -F"\t" -v trait="$trait_name" -v ss_tol="$ss_tol" '
        #    BEGIN { OFS="\t"; ss_tol += 0 }
        #    NR==1 { for(i=1;i<=NF;i++) if($i==ss_col) ss_idx=i; next }
        #    $ss_idx!="" && $ss_idx ~ /^[0-9.]+$/ { vals[++n]=$ss_idx }
        #    END {
        #        if(n > 0) {
        #            asort(vals)
        #            if(n%2) {med=vals[int(n/2)+1]}
        #            else {med=(vals[int(n/2)]+vals[int(n/2)+1])/2}
        #            low  = med * (1 - ss_tol)
        #            high = med * (1 + ss_tol)
        #            print trait, low, med, high
        #        } else {
        #            print trait, "NA", "NA", "NA"
        #        }
        #    }
        #' >> "$trait_summary_table"

        # extract numeric SS values, sort numerically (disk-backed), get count
        tmp="$workdir/${trait_name}.ssvals.sorted"

	make_ss_data > "$tmp"

        if [[ -s "$tmp" ]]; then
            n=$(wc -l < "$tmp")
            med=$(awk -v n="$n" '
                NR==int((n+1)/2){a=$1}
                NR==int((n+2)/2){b=$1}
                END{print (n%2)?a:(a+b)/2}
            ' "$tmp")

           low=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1-t)}')
           high=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1+t)}')

           printf "%s\t%s\t%s\t%s\n" "$trait_name" "$low" "$med" "$high" >> "$trait_summary_table"
        else
           printf "%s\tNA\tNA\tNA\n" "$trait_name" >> "$trait_summary_table"
        fi

        rm -f "$tmp"

        echo "created trait ss table"
        memsnap
    fi
done

echo "finished trait loop"
memsnap

kill "$sampler" 2>/dev/null || true
wait "$sampler" 2>/dev/null || true
