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
workdir="1_comb_form_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$workdir"

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
#col_pos=$(get_col "pos")
col_chrom=$(get_col "chrom")
col_A1=$(get_col "A1")
col_A2=$(get_col "A2")
col_beta_hat=$(get_col "beta_hat")
col_se=$(get_col "se")
#col_pval=$(get_col "p_value")
col_af=$(get_col "af")
col_sample_size=$(get_col "sample_size")
#col_pub_sample_size=$(get_col "pub_sample_size")
#col_effect_is_or=$(get_col "effect_is_or")
col_name=$(get_col "name")

# loop over gwas files and format - - -

trait_summary_table="$workdir/trait_sample_stats.tsv"
# here common means SNPs in common between all traits.  this will confuse people bc of af so I need to change
common_snps_and_maf="$workdir/common_snps_and_maf.tsv"
common_snps_and_maf_tmp="$workdir/common_snps_and_maf.tmp"

# Write header to summary table (only once)
echo -e "trait\tss_low\tss_median\tss_high" > "$trait_summary_table" 

for ((i=2; i<=num_traits+1; i++)); do
#for ((i=2; i<=2; i++)); do
#for ((i=2; i<=3; i++)); do
    # Read info fields for trait i (awk column indexing, adjust as needed for TSV)
    f=$(awk -F, -v row="$i" -v col="$col_raw_data" 'NR==row {print $col}' "$gwas_info_file")
    snp=$(awk -F, -v row="$i" -v col="$col_snp" 'NR==row {print $col}' "$gwas_info_file")
#    pos=$(awk -F, -v row="$i" -v col="$col_pos" 'NR==row {print $col}' "$gwas_info_file")
    chrn=$(awk -F, -v row="$i" -v col="$col_chrom" 'NR==row {print $col}' "$gwas_info_file")
    A1=$(awk -F, -v row="$i" -v col="$col_A1" 'NR==row {print $col}' "$gwas_info_file")
    A2=$(awk -F, -v row="$i" -v col="$col_A2" 'NR==row {print $col}' "$gwas_info_file")
    beta_hat=$(awk -F, -v row="$i" -v col="$col_beta_hat" 'NR==row {print $col}' "$gwas_info_file")
    se=$(awk -F, -v row="$i" -v col="$col_se" 'NR==row {print $col}' "$gwas_info_file")
#    pval=$(awk -F, -v row="$i" -v col="$col_pval" 'NR==row {print $col}' "$gwas_info_file")
    af=$(awk -F, -v row="$i" -v col="$col_af" 'NR==row {print $col}' "$gwas_info_file")
    sample_size=$(awk -F, -v row="$i" -v col="$col_sample_size" 'NR==row {print $col}' "$gwas_info_file")
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
        upstream() {
 		 zcat "$f" \
  			| awk -F"$delimiter" -v col="$chr_col" -v chrom_val="$chrom" -v OFS="\t" 'NR==1 || $col == chrom_val' \
  			| awk -F"\t" -v OFS="\t" \
      				-v snp_name="$snp" -v A1_name="$A1" -v A2_name="$A2" \
				-v beta_name="$beta_hat" -v se_name="$se" \
      				-v ss_name="$sample_size" -v af_name="$af" \
      				-f remove_invalid_variants.awk \
  			| awk -F"\t" -v OFS="\t" -v beta_name="$beta_hat" -v se_name="$se" -f add_zscore.awk
	}


        upstream | wc -l  # CHECK IF NOT TRUNCATED

        #upstream 2>/dev/null | head

	# this fifo thing is named pipe to help us not run upstream twice (as you would to just use paste).  I'm not sure it would work well for ppl locally.  but reading 2x is slow
        # keep the fifo for now
	#make_filt_data() {
        #    local tmpdir fifo
        #    tmpdir=$(mktemp -d)
        #    fifo="$tmpdir/pvals.fifo"
        #    mkfifo "$fifo"
            
	#    cleanup() { rm -rf "$tmpdir"; }
        #    trap cleanup EXIT INT TERM

        #    upstream \
        #    | tee >(  # side branch computes P and writes to fifo
        #        awk -F$'\t' '
        #            NR==1 {
        #                for (i=1; i<=NF; i++) if ($i=="Z") { z=i; break }
        #                if (!z) { print "ERROR: no Z column in header" > "/dev/stderr"; exit 2 }
        #                next
        #            }
        #            { print $z }
        #        ' \
        #    | Rscript -e '
        #        z_chr <- scan("stdin", what=character(), quiet=TRUE)
        #        z_num <- suppressWarnings(as.numeric(z_chr))
        #        p <- 2*pnorm(-abs(z_num))
        #        out <- ifelse(is.na(z_num), "NA", format(p, scientific=TRUE, digits=16))
        #        cat("P\n"); cat(out, sep="\n")
        #    ' \
        #    > "$fifo"
        #    ) \
        #    | paste -d $'\t' - "$fifo"

            # cleanup will run via trap
       # }
	
	
	make_filt_data() {
  	    paste -d $'\t' \
    	        <( upstream ) \
    		<(
      		    upstream \
      		        | awk -F$'\t' '
          		    NR==1 {
            			for (i=1; i<=NF; i++) if ($i=="Z") { z=i; break }
            			if (!z) { print "ERROR: no Z column in header" > "/dev/stderr"; exit 2 }
            			next
          			}
          		    { print $z }
        		    ' \
      			| Rscript -e '
          		    z_chr <- scan("stdin", what=character(), quiet=TRUE)
          		    z_num <- suppressWarnings(as.numeric(z_chr))
          		    p <- 2*pnorm(-abs(z_num))
          		    out <- ifelse(is.na(z_num), "NA", format(p, scientific=TRUE, digits=16))
          		    cat("p\n"); cat(out, sep="\n")
        		'
   		 )
	}  
	    
    	echo "added p-values to filtered data from R"


	#make_filt_data | wc -l  # CHECK IF NOT TRUNCATED
        #make_filt_data 2>/dev/null | head

	# 1. Sample size statistics for trait summary table
	make_filt_data | awk -F"\t" -v ss_col="$sample_size" -v trait="$trait_name" '
            BEGIN { OFS="\t" }
	    NR==1 { for(i=1;i<=NF;i++) if($i==ss_col) ss_idx=i; next }
            $ss_idx!="" && $ss_idx ~ /^[0-9.]+$/ { vals[++n]=$ss_idx }
            END {
                if(n > 0) {
                    asort(vals)
                    if(n%2) {med=vals[int(n/2)+1]}
                    else {med=(vals[int(n/2)]+vals[int(n/2)+1])/2}
                    low=med*0.9
                    high=med*1.1
                    print trait, low, med, high
                } else {
                    print trait, "NA", "NA", "NA"
                }
            }
        ' >> "$trait_summary_table"

	echo "created trait ss table"

	read trait low med high < <(awk -F"\t" -v trait="$trait_name" '$1==trait {print $1, $2, $3, $4}' "$trait_summary_table" | tail -n1)

	# 2. Common SNPs and min MAF intersection/update
        if ((i == 2)); then
            make_filt_data | \
	    awk -F"\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" -v p_col="p" -v low="$low" -v high="$high" '
	        BEGIN { OFS="\t" }
                NR==1 {
                    for (i=1; i<=NF; i++) {
                        if ($i == snp_col) snp_idx = i
                        if ($i == af_col) af_idx = i
                        if ($i == ss_col) ss_idx = i
			#if ($i == z_col) z_idx = i
			if ($i == p_col) p_idx = i
                    }
                    next
                }
                $snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
                    maf = ($af_idx < 1 - $af_idx) ? $af_idx : 1 - $af_idx
                    ss_val = $ss_idx + 0
                    ss_flag = (ss_val >= low && ss_val <= high) ? 1 : 0
		    #z = (z_idx && $(z_idx)!="") ? $(z_idx) : "NA"
		    p = (p_idx && $(p_idx)!="") ? $(p_idx) : "NA"
                    print $snp_idx, maf, p, ss_flag
                }
                ' > "$common_snps_and_maf"
        else
	    make_filt_data | \
            awk -F"\t" -v snp_col="$snp" \
                -v af_col="$af" \
                -v ss_col="$sample_size" \
                -v p_col="p" \
		-v low="$low" \
                -v high="$high" \
                'BEGIN { OFS="\t" }
                # First file: filtered_data, build lookup tables
                NR==FNR {
                    if(FNR==1) {
                        for(i=1;i<=NF;i++) {
                            if($i==snp_col) snp_idx=i
                            if($i==af_col) af_idx=i
			    if($i == p_col) p_idx = i
                            if($i==ss_col) ss_idx=i
                        }
                        #print "Indices:", snp_idx, af_idx, ss_idx;
                        next
		    } else {
        	        # Only store filtered_data in lookup
                        af_val = $af_idx + 0
			p_val = $p_idx + 0
                        ss_val = $ss_idx + 0
                        snp_val = $snp_idx
                        filtered_af[snp_val] = af_val
                        filtered_p[snp_val] = p_val
			filtered_ss[snp_val] = ss_val
                        #print "Raw snp:", $snp_idx, "Raw af:", $af_idx, "Raw ss:", $ss_idx
		        #print "Storing:", snp_val, "af:", af_val, "ss:", ss_val
	            }
		}
                # Second file: common_snps_and_maf
                NR!=FNR {
                    snp = $1
                    prior_maf = $2 + 0
		    prior_p = $3 + 0
                    flag = $4
                    # Only keep SNPs present in filtered_data
                    if (snp in filtered_af && snp in filtered_p && snp in filtered_ss) {
                        # Choose minimum MAF
                        maf1 = filtered_af[snp]
                        maf2 = 1 - filtered_af[snp]
                        min_maf = prior_maf
	       		if(maf1 < min_maf) min_maf = maf1
                        if(maf2 < min_maf) min_maf = maf2
			# Choose minimum p
                        p1 = filtered_p[snp]
                        min_p = prior_p
                        if(p1 < min_p) min_p = p1
			# Evaluate sample size flag
                        ss_val = filtered_ss[snp]
                        new_flag = flag
                        if(flag==1 && (ss_val < low || ss_val > high)) new_flag = 0
                        print snp, min_maf, min_p, new_flag
                    }
               }
               ' - "$common_snps_and_maf" > "$common_snps_and_maf_tmp" &&
               mv "$common_snps_and_maf_tmp" "$common_snps_and_maf"
		
        echo "created comm snps table"

        fi
    fi
done

# ---
awk -F"\t" 'BEGIN { OFS="\t"; print "snp", "min_maf", "min_p", "in_ss_range_each_trait", "above_min_maf_thresh" }
{
    af_flag = ($2 >= 0.05) ? 1 : 0
    print $1, $2, $3, $4, af_flag
}' "$common_snps_and_maf" > "$common_snps_and_maf_tmp" &&
mv "$common_snps_and_maf_tmp" "$common_snps_and_maf"

echo "wrote out snp table with maf and ss flags"

# write out a file where only the snps that pass ss and maf filters stay
final_pass_snps="$workdir/snps_pass_all_filts.txt"
awk -F"\t" 'NR>1 && $4==1 && $5==1 {print $1}' "$common_snps_and_maf" > "$final_pass_snps"

echo "wrote out passing snps:  found in all traits, pass ss filter, pass maf filter"

echo "finished formatting"
