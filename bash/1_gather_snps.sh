#!/bin/bash

set -x

set -euo pipefail

# make all awks, joins, and sorts consistent
export LC_ALL=C

# temp workspace for generated files
workdir="1_comb_form_$(date +%Y%m%d_%H%M%S)"
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
  done ) > "$workdir/1_ps_rss.log" &
sampler=$!

# read input arguments - - -

#bash bash/1_combine_and_format.sh {wildcards.chrom} {input.gwas_info} {input.ref_bim} {params.af_thresh} {params.sample_size_tol} {output.out}

chrom="$1"                # e.g., from Snakemake wildcards
gwas_info_file="$2"       # path to info file
trait_summary_table="$3"  # table of ss medians for each trait
af_thresh="$4"            # allele freq threshold
ss_tol="$5"      # sample size tolerance
out="$6"                  # output file

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

#trait_summary_table="$workdir/trait_sample_stats.tsv"
# here shared means SNPs in shared between all traits. 
shared_snps_and_maf="$workdir/shared_snps_and_maf.tsv"
shared_snps_and_maf_tmp="$workdir/shared_snps_and_maf.tmp"

# Write header to summary table (only once)
#echo -e "trait\tss_low\tss_median\tss_high" > "$trait_summary_table" 

for ((i=2; i<=num_traits+1; i++)); do
#for ((i=3; i<=3; i++)); do
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

	make_filt_data() {
            zcat "$f" |
                awk -F"$delimiter" -v chr_header="$chrn" -v chrom_val="$chrom" -v OFS="\t" '
                    NR==1 {
                        for (i=1; i<=NF; i++) if ($i == chr_header) chr_idx = i
                        if (!chr_idx) { print "ERROR: chrom header not found: " chr_header > "/dev/stderr"; exit 1 }
                        print
                        next
                    }
                    $chr_idx == chrom_val
                    ' \
		| awk -F"\t" -v OFS="\t" \
      			-v snp_name="$snp" -v A1_name="$A1" -v A2_name="$A2" \
			-v beta_name="$beta_hat" -v se_name="$se" \
      			-v ss_name="$sample_size" -v af_name="$af" \
      			-f remove_invalid_variants.awk \
		| { read -r header; printf '%s\n' "$header"; sort -T "$workdir" -S 100M -t $'\t' -k1,1 -u; } \
  		| awk -F"\t" -v OFS="\t" -v beta_name="$beta_hat" -v se_name="$se" -f add_zscore.awk
	}
        # add the copying over of pub ss if real ss not avail
   
	#make_filt_data > tmp.txt
	#echo "maxrss after using make_filt_data"
        #memsnap

	#make_filt_data | wc -l  # CHECK IF NOT TRUNCATED
        #make_filt_data 2>/dev/null | head

	# 1. Sample size statistics for trait summary table
	#read trait low med high < <(awk -F"\t" -v trait="$trait_name" '$1==trait {print $1, $2, $3, $4}' "$trait_summary_table" | tail -n1)
        read -r trait med < <(awk -F"\t" -v trait="$trait_name" '$1==trait {print $1, $2}' "$trait_summary_table")
	low=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1-t)}')
        high=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1+t)}')

        echo $low
	echo $high

        #low=1000
	#high=10000

	# 2. Common SNPs and min MAF intersection/update
        if ((i == 2)); then
            make_filt_data | \
	    awk -F"\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" -v low="$low" -v high="$high" '
	        BEGIN { OFS="\t" }
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        if ($i == snp_col) snp_idx = i
                        if ($i == af_col) af_idx = i
                        if ($i == ss_col) ss_idx = i
			if ($i == z_col) z_idx = i
                    }
                    next
                }
                $snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
                    maf = ($af_idx < 1 - $af_idx) ? $af_idx : 1 - $af_idx
                    ss_val = $ss_idx + 0
                    ss_flag = (ss_val >= low && ss_val <= high) ? 1 : 0
		    # absolute the z
		    z = (z_idx && $(z_idx)!="") ? ($(z_idx)+0) : "NA"
		    if (z != "NA" && z < 0) z = -z
                    print $snp_idx, maf, z, ss_flag
                }
            ' > "$shared_snps_and_maf"
	    echo "after making first shared snp table"
	    memsnap
        else
            #trait_keyed="$(mktemp)"
            #trait_sorted="$(mktemp)"

            # process sub extract SNP, AF, SS, |Z| from this trait (header-aware), headerless output
                        		
            # INNER JOIN shared (SNP prior_maf prior_z ss_flag) with trait (SNP af ss zabs)
            # Only SNPs present in BOTH will be output => intersection across traits

	    join -t $'\t' -1 1 -2 1 \
                -o 1.1,1.2,1.3,1.4,2.2,2.3,2.4 \
                "$shared_snps_and_maf" \
                <(
                  make_filt_data |
                    awk -F"\t" -v OFS="\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" '
                        NR==1{
                            for(i=1;i<=NF;i++){
                                if($i==snp_col) snp_idx=i
                                if($i==af_col)  af_idx=i
                                if($i==ss_col)  ss_idx=i
                                if($i==z_col)   z_idx=i
                            }
                            next
                        }
                        $snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
                        z = (z_idx && $(z_idx)!="") ? ($(z_idx)+0) : "NA"
                        if(z!="NA" && z<0) z=-z
                        print $snp_idx, ($af_idx+0), ($ss_idx+0), z
                        }
                    '
                 ) |
            awk -F"\t" -v OFS="\t" -v low="$low" -v high="$high" '
                {
                    snp=$1
                    prior_maf=$2+0
                    prior_z=$3
                    ss_flag=$4+0
                    af=$5+0
                    ss=$6+0
                    z=$7

                    # update min_maf using af and 1-af
                    maf2=1-af
                    min_maf=prior_maf
                    if(af < min_maf)   min_maf=af
                    if(maf2 < min_maf) min_maf=maf2

                    # update max_abs_z (treat NA as missing)
                    max_abs_z=prior_z
                    if(max_abs_z!="NA"){
                        max_abs_z+=0
                        if(max_abs_z<0) max_abs_z=-max_abs_z
                    }
                    if(z!="NA"){
                        z+=0
                        if(z<0) z=-z
                        if(max_abs_z=="NA" || z>max_abs_z) max_abs_z=z
                    }

                    new_ss_flag=ss_flag
                    if(ss_flag==1 && (ss < low || ss > high)) new_ss_flag=0

                    print snp, min_maf, max_abs_z, new_ss_flag
                }
            ' > "$shared_snps_and_maf_tmp"

            mv "$shared_snps_and_maf_tmp" "$shared_snps_and_maf"

	    echo "updated shared snps table"
            memsnap

            #rm -f "$trait_keyed" "$trait_sorted"
        fi
    fi
done

echo "finished trait loop"
memsnap

# ---
awk -F"\t" -v af_thresh="$af_thresh" '
BEGIN { OFS="\t"; print "snp", "min_maf", "max_abs_z", "in_ss_range_each_trait", "above_min_maf_thresh" }
{
    af_flag = ($2 >= af_thresh) ? 1 : 0
    print $1, $2, $3, $4, af_flag
}' "$shared_snps_and_maf" > "$shared_snps_and_maf_tmp" &&
mv "$shared_snps_and_maf_tmp" "$shared_snps_and_maf"

echo "wrote out snp table with maf and ss flags"
memsnap

# write out a file where only the snps that pass ss and maf filters stay
final_pass_snps="$workdir/$out"
awk -F"\t" 'NR==1 {print $1"\t"$3; next} $4==1 && $5==1 {print $1"\t"$3}' "$shared_snps_and_maf" > "$final_pass_snps"
awk -F"\t" 'NR>1 && $4==1 && $5==1 {print $1}' "$shared_snps_and_maf" > "$final_pass_snps"

echo "wrote out passing snps:  found in all traits, pass ss filter, pass maf filter"
memsnap

echo "finished formatting"
memsnap

kill "$sampler" 2>/dev/null || true
wait "$sampler" 2>/dev/null || true
