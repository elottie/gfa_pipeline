#!/bin/bash

#set -x

set -euo pipefail

memsnap() {
  # prints MaxRSS and AveRSS for the current job step
  sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxRSS,AveRSS -n 2>/dev/null \
    | awk '{print strftime("[%F %T]"), $0}'
}

stat -fc %T /sys/fs/cgroup

# --- cgroup memory introspection (works on most cgroup v2 systems) ---
cgmem_snap () {
  local tag="${1:-snap}"
  local cg="/sys/fs/cgroup"

  echo "=== [$tag] $(date) ==="

  # SLURM context (helps interpret sacct output)
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-} SLURM_STEP_ID=${SLURM_STEP_ID:-} SLURM_PROCID=${SLURM_PROCID:-}"

  if [[ -r "$cg/memory.current" ]]; then
    echo "cgroup: $(cat $cg/cgroup.controllers 2>/dev/null | tr -s ' ' | head -c 200) ..."
    echo "memory.current: $(cat $cg/memory.current)"
    echo "memory.max:     $(cat $cg/memory.max)"
  else
    echo "No $cg/memory.current (maybe cgroup v1 or different mount); skipping cgroup v2 stats."
  fi

  if [[ -r "$cg/memory.stat" ]]; then
    # Key fields:
    #  anon = heap/stack (your process memory)
    #  file = page cache charged to cgroup
    #  slab = kernel objects (dentries/inodes etc.)
    awk '
      BEGIN { want["anon"]=1; want["file"]=1; want["slab"]=1;
              want["active_file"]=1; want["inactive_file"]=1;
              want["active_anon"]=1; want["inactive_anon"]=1 }
      $1 in want { printf "%-14s %s\n", $1":", $2 }
    ' "$cg/memory.stat"
  fi

  echo
}

# Optional: per-process RSS of the current bash (won't capture file cache)
procrss_snap () {
  local tag="${1:-rss}"
  echo "=== [$tag] bash RSS ==="
  awk '/VmRSS|VmHWM/ {print}' /proc/$$/status
  echo
}

( while :; do
    date
    ps -u "$USER" -o pid,ppid,rss,cmd --sort=-rss | head -20
    echo "----"
    sleep 2
  done ) > ps_rss.log &
sampler=$!

# read input arguments - - -

#bash bash/1_combine_and_format.sh {wildcards.chrom} {input.gwas_info} {input.ref_bim} {params.af_thresh} {params.sample_size_tol} {output.out}

chrom="$1"                # e.g., from Snakemake wildcards
gwas_info_file="$2"       # path to info file
af_thresh="$3"            # allele freq threshold
ss_tol="$4"      # sample size tolerance
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

echo "after extracting col names"
memsnap

# loop over gwas files and format - - -

trait_summary_table="$workdir/trait_sample_stats.tsv"
# here shared means SNPs in shared between all traits. 
shared_snps_and_maf="$workdir/shared_snps_and_maf.tsv"
shared_snps_and_maf_tmp="$workdir/shared_snps_and_maf.tmp"

# Write header to summary table (only once)
echo -e "trait\tss_low\tss_median\tss_high" > "$trait_summary_table" 

#for ((i=2; i<=num_traits+1; i++)); do
#for ((i=3; i<=3; i++)); do
for ((i=2; i<=3; i++)); do
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

    delimiter=$(get_file_delimiter "$f")

    #trait_out="$workdir/${trait_name}.final.tsv"
    if [[ "$f" == *.vcf.gz || "$f" == *.vcf.bgz ]]; then
        echo "Calling format_ieu_chrom (external): $f $chrom $af_thresh"
    #    format_ieu_chrom "$f" "$chrom" "$af_thresh" > "$trait_out"
    else
        echo "maxrss before make_filt_data"
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
   
	cgmem_snap "iter=$i BEFORE make_filt_data"
        procrss_snap "iter=$i BEFORE make_filt_data"
	
	make_filt_data > tmp.txt
	
        cgmem_snap "iter=$i AFTER make_filt_data"
        procrss_snap "iter=$i AFTER make_filt_data"

	echo "maxrss after make_filt_data"
        memsnap

#	echo "----- iter $i -----"
#	sstat -j "$SLURM_JOB_ID" --format=JobID,MaxRSS,AveRSS -n

	# show cgroup v1 cache vs rss if available
	#cgbase=/sys/fs/cgroup/memory
	#jobcg=$(find "$cgbase" -maxdepth 6 -type d -name "job_${SLURM_JOB_ID}*" 2>/dev/null | head -n 1)
	#echo "jobcg=$jobcg"
	#if [[ -n "$jobcg" && -r "$jobcg/memory.stat" ]]; then
  #		egrep '^(rss|cache|mapped_file|inactive_file|active_file|inactive_anon|active_anon) ' "$jobcg/memory.stat"
#	fi

	#printf "A1=%s A2=%s f=%s chrom=%s delimiter=%q chrn=%s snp=%s beta_hat=%s se=%s sample_size=%s af=%s\n" \
        # "$A1" "$A2" "$f" "$chrom" "$delimiter" "$chrn" "$snp" "$beta_hat" "$se" "$sample_size" "$af"

	#export -f make_filt_data
        #export f delimiter chrn chrom snp A1 A2 beta_hat se sample_size af  # and any others it references
        #/usr/bin/time -v bash -lc 'make_filt_data > /dev/null'

            #/usr/bin/time -v env LC_ALL=C sort -S 200M -t $'\t' -k1,1 "$trait_keyed" > /dev/null
	
	#n=$(make_filt_data | wc -l)
#	echo "trait=$trait_name make_filt_data_lines=$n" >&2

	#/usr/bin/time -v bash -c 'make_filt_data > /dev/null'

	#ps -eo pid,ppid,rss,cmd --sort=-rss | head -20

	#make_filt_data | wc -l  # CHECK IF NOT TRUNCATED
        #make_filt_data 2>/dev/null | head

	# 1. Sample size statistics for trait summary table
	#make_filt_data | awk -F"\t" -v ss_col="$sample_size" -v trait="$trait_name" -v ss_tol="$ss_tol" '
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
	#tmp="$workdir/${trait_name}.ssvals"
	#make_filt_data | awk -F"\t" -v ss_col="$sample_size" '
 	#     NR==1 {for(i=1;i<=NF;i++) if($i==ss_col) ss=i; next}
        #     $ss!="" && $ss ~ /^[0-9.]+$/ {print $ss}
        #     ' > "$tmp"

        #n=$(wc -l < "$tmp")
        #if (( n > 0 )); then
        #    med=$(sort -n "$tmp" | awk -v n="$n" '
        #        NR==int((n+1)/2){a=$1}
        #        NR==int((n+2)/2){b=$1}
        #        END{print (n%2)?a:(a+b)/2}
        #    ')
        #    low=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1-t)}')
        #    high=$(awk -v m="$med" -v t="$ss_tol" 'BEGIN{print m*(1+t)}')
        #    printf "%s\t%s\t%s\t%s\n" "$trait_name" "$low" "$med" "$high" >> "$trait_summary_table"
        #else
        #    printf "%s\tNA\tNA\tNA\n" "$trait_name" >> "$trait_summary_table"
        #fi
        #rm -f "$tmp"

	#echo "created trait ss table"
        #memsnap

	#read trait low med high < <(awk -F"\t" -v trait="$trait_name" '$1==trait {print $1, $2, $3, $4}' "$trait_summary_table" | tail -n1)
        low=1000
	high=10000

	# 2. Common SNPs and min MAF intersection/update
        #if ((i == 2)); then
        #    make_filt_data | \
	#    awk -F"\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" -v low="$low" -v high="$high" '
	#        BEGIN { OFS="\t" }
        #        NR==1 {
        #            for (i=1; i<=NF; i++) {
        #                if ($i == snp_col) snp_idx = i
        #                if ($i == af_col) af_idx = i
        #                if ($i == ss_col) ss_idx = i
#			if ($i == z_col) z_idx = i
 #                   }
 #                   next
 #               }
 #               $snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
 #                   maf = ($af_idx < 1 - $af_idx) ? $af_idx : 1 - $af_idx
 #                   ss_val = $ss_idx + 0
 #                   ss_flag = (ss_val >= low && ss_val <= high) ? 1 : 0
#		    # absolute the z
#		    z = (z_idx && $(z_idx)!="") ? ($(z_idx)+0) : "NA"
#		    if (z != "NA" && z < 0) z = -z
#                    print $snp_idx, maf, z, ss_flag
#                }
#            ' | LC_ALL=C sort -S 200M -t $'\t' -k1,1 > "$shared_snps_and_maf"
#	    echo "after making first shared snp table"
#	    memsnap
#	    export -f make_filt_data
#	    export snp af sample_size low high shared_snps_and_maf f delimiter chrn chrom snp A1 A2 beta_hat se sample_size af # plus any vars make_filt_data uses

#	/usr/bin/time -v bash -lc '
#  	make_filt_data |
#  	awk -F"\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" -v low="$low" -v high="$high" '"'"'
#    	BEGIN { OFS="\t" }
#    	NR==1 {
#      		for (i=1; i<=NF; i++) {
#       		 if ($i == snp_col) snp_idx = i
#        	if ($i == af_col) af_idx = i
#        	if ($i == ss_col) ss_idx = i
#        	if ($i == z_col) z_idx = i
#      	}
#      	next
#    	}
#    	$snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
#      	maf = ($af_idx < 1 - $af_idx) ? $af_idx : 1 - $af_idx
#      	ss_val = $ss_idx + 0
#      	ss_flag = (ss_val >= low && ss_val <= high) ? 1 : 0
#      	z = (z_idx && $(z_idx)!="") ? ($(z_idx)+0) : "NA"
#      	if (z != "NA" && z < 0) z = -z
#      	print $snp_idx, maf, z, ss_flag
#    	}
#  	'"'"' |
#  	LC_ALL=C sort -S 200M -t $'"'"'\t'"'"' -k1,1 > "$shared_snps_and_maf"
#	'

#       else
	#    make_filt_data | \
        #    awk -F"\t" -v snp_col="$snp" \
        #        -v af_col="$af" \
        #        -v ss_col="$sample_size" \
        #        -v z_col="Z" \
	#	-v low="$low" \
        #        -v high="$high" \
        #        'BEGIN { OFS="\t" }
        #        # First file: filtered_data, build lookup tables
        #        NR==FNR {
        #            if(FNR==1) {
        #                for(i=1;i<=NF;i++) {
        #                    if($i==snp_col) snp_idx=i
        #                    if($i==af_col) af_idx=i
        #                    if($i==ss_col) ss_idx=i
	# 	            if($i == z_col) z_idx = i
        #                }
        #                #print "Indices:", snp_idx, af_idx, ss_idx;
        #                next
	#		    } else {
        # 	        # Only store filtered_data in lookup
        #                snp_val = $snp_idx
	#			af_val = $af_idx + 0
        #                ss_val = $ss_idx + 0
	#		z_val = $z_idx + 0
        #                filtered_af[snp_val] = af_val
	#		filtered_ss[snp_val] = ss_val
#			filtered_z[snp_val] = z_val
#                        #print "Raw snp:", $snp_idx, "Raw af:", $af_idx, "Raw ss:", $ss_idx
##		        #print "Storing:", snp_val, "af:", af_val, "ss:", ss_val
#	            }
#		}
#                # Second file: shared_snps_and_maf
#                NR!=FNR {
#                    snp = $1
#                    prior_maf = $2 + 0
#		    prior_z = $3 + 0
#                    ss_flag = $4
#                    # Only keep SNPs present in filtered_data
#                    if (snp in filtered_af && snp in filtered_ss && snp in filtered_z) {
#                        
#			# Choose minimum MAF
#                        maf1 = filtered_af[snp]
#                        maf2 = 1 - filtered_af[snp]
#                        min_maf = prior_maf
#	       		if(maf1 < min_maf) min_maf = maf1
#                        if(maf2 < min_maf) min_maf = maf2
#			
#			# Choose absolute maximum Z
#                        z1 = filtered_z[snp]
#			# absolute the prior_z
#			max_abs_z = prior_z
#			if (max_abs_z < 0) max_abs_z = -max_abs_z
#			# absolute the current z
#			abs_z1 = z1
#			if (abs_z1 < 0) abs_z1 = -abs_z1
#			# choose the greater absolute z
#			if (abs_z1 > max_abs_z) max_abs_z = abs_z1
#			
#			# Evaluate sample size flag
#                        ss_val = filtered_ss[snp]
#                        new_ss_flag = ss_flag
#                        if(ss_flag==1 && (ss_val < low || ss_val > high)) new_ss_flag = 0
#                       
#			print snp, min_maf, max_abs_z, new_ss_flag
#                    }
#               }
#               ' - "$shared_snps_and_maf" > "$shared_snps_and_maf_tmp" &&
#               mv "$shared_snps_and_maf_tmp" "$shared_snps_and_maf"

	#else
 #           trait_keyed="$(mktemp)"
#            trait_sorted="$(mktemp)"

            # Extract SNP, AF, SS, |Z| from this trait (header-aware), headerless output
#            make_filt_data |
#                awk -F"\t" -v OFS="\t" -v snp_col="$snp" -v af_col="$af" -v ss_col="$sample_size" -v z_col="Z" '
#                    NR==1{
#                        for(i=1;i<=NF;i++){
#                            if($i==snp_col) snp_idx=i
#                            if($i==af_col)  af_idx=i
#                            if($i==ss_col)  ss_idx=i
#                            if($i==z_col)    z_idx=i
#                        }
#                        next
#                    }
#                    $snp_idx!="" && $af_idx!="" && $af_idx ~ /^[0-9.]+$/ && $af_idx>=0 && $af_idx<=1 && $ss_idx ~ /^[0-9.]+$/{
#                        z = (z_idx && $(z_idx)!="") ? ($(z_idx)+0) : "NA"
#                        if(z!="NA" && z<0) z=-z
#                        print $snp_idx, ($af_idx+0), ($ss_idx+0), z
#                    }
#               ' > "$trait_keyed"
#            LC_ALL=C sort -S 200M -t $'\t' -k1,1 "$trait_keyed" > "$trait_sorted"
		
#	    export -f make_filt_data
#	    export f delimiter chrn chrom snp A1 A2 beta_hat se sample_size af prefix  # and any others it references
#	    /usr/bin/time -v bash -lc 'make_filt_data > /dev/null'
#		
#	    /usr/bin/time -v env LC_ALL=C sort -S 200M -t $'\t' -k1,1 "$trait_keyed" > /dev/null

            # INNER JOIN shared (SNP prior_maf prior_z ss_flag) with trait (SNP af ss zabs)
            # Only SNPs present in BOTH will be output => intersection across traits
#            join -t $'\t' -1 1 -2 1 \
#                -o 1.1,1.2,1.3,1.4,2.2,2.3,2.4 \
#                "$shared_snps_and_maf" "$trait_sorted" |
#            awk -F"\t" -v OFS="\t" -v low="$low" -v high="$high" '
#                {
#                    snp=$1
#                    prior_maf=$2+0
#                    prior_z=$3
#                    ss_flag=$4+0
#                    af=$5+0
#                    ss=$6+0
#                    z=$7

                    # update min_maf using af and 1-af
#                    maf2=1-af
#                    min_maf=prior_maf
#                    if(af < min_maf)   min_maf=af
#                    if(maf2 < min_maf) min_maf=maf2

                    # update max_abs_z (treat NA as missing)
#                    max_abs_z=prior_z
#                    if(max_abs_z!="NA"){
#                        max_abs_z+=0
#                        if(max_abs_z<0) max_abs_z=-max_abs_z
#                    }
#                    if(z!="NA"){
#                        z+=0
#                        if(z<0) z=-z
#                        if(max_abs_z=="NA" || z>max_abs_z) max_abs_z=z
#                    }

#                    new_ss_flag=ss_flag
#                    if(ss_flag==1 && (ss < low || ss > high)) new_ss_flag=0

#                    print snp, min_maf, max_abs_z, new_ss_flag
#                }
#            ' > "$shared_snps_and_maf_tmp"

#            mv "$shared_snps_and_maf_tmp" "$shared_snps_and_maf"

	    #ps -eo pid,ppid,rss,cmd --sort=-rss | head -20

#	    echo "updated shared snps table"
#            memsnap

#            rm -f "$trait_keyed" "$trait_sorted"
#        fi
    #echo "END ITER $i: $(date)"
    #jobs -l || true
    #ps -u "$USER" -o pid,ppid,stat,cmd --sort=ppid | awk '$4 ~ /zcat|awk|sort|gzip/ {print}'
    #wait
    #kill %1 2>/dev/null
   
    echo "END ITER 1 $(date)"
    #jobs -l || true
    #ps -u "$USER" -o pid,ppid,stat,rss,cmd --sort=-rss | head -25
    #wait
    #echo "AFTER WAIT $(date)"
    #ps -u "$USER" -o pid,ppid,stat,rss,cmd --sort=-rss | head -10


    fi
done

echo "finished trait loop"
memsnap

# ---
#awk -F"\t" -v af_thresh="$af_thresh" '
#BEGIN { OFS="\t"; print "snp", "min_maf", "max_abs_z", "in_ss_range_each_trait", "above_min_maf_thresh" }
#{
#    af_flag = ($2 >= af_thresh) ? 1 : 0
#    print $1, $2, $3, $4, af_flag
#}' "$shared_snps_and_maf" > "$shared_snps_and_maf_tmp" &&
#mv "$shared_snps_and_maf_tmp" "$shared_snps_and_maf"

#echo "wrote out snp table with maf and ss flags"
#memsnap

# write out a file where only the snps that pass ss and maf filters stay
#final_pass_snps="$workdir/$out"
#awk -F"\t" 'NR==1 {print $1"\t"$3; next} $4==1 && $5==1 {print $1"\t"$3}' "$shared_snps_and_maf" > "$final_pass_snps"
#awk -F"\t" 'NR>1 && $4==1 && $5==1 {print $1}' "$shared_snps_and_maf" > "$final_pass_snps"

#echo "wrote out passing snps:  found in all traits, pass ss filter, pass maf filter"
#memsnap

#echo "finished formatting"
#memsnap

kill "$sampler" 2>/dev/null || true
wait "$sampler" 2>/dev/null || true
