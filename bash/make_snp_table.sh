#!/bin/bash

gwas_info_file="../C100001554_And_Friends_3Metabolites.csv"

read highest_min lowest_max < <(
    awk -F',' 'NR==1 {
        for(i=1;i<=NF;i++) {
            if($i=="raw_data_path") f_col=i;
            if($i=="sample_size") s_col=i;
            if($i=="chrom") c_col=i;
        }
        next
    }
    {print $f_col "\t" $s_col "\t" $c_col}
    ' "$gwas_info_file" \
    | while IFS=$'\t' read file sample_col chrom_col; do
        [[ "$file" == *.gz ]] && cat_cmd="zcat" || cat_cmd="cat"

        # Get column numbers for sample size and chromosome
        sample_num=$($cat_cmd "$file" | awk -F'\t' -v col="$sample_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')
        chr_num=$($cat_cmd "$file" | awk -F'\t' -v col="$chrom_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')

        # Compute median for chromosome 19 sample sizes
        $cat_cmd "$file" | awk -F'\t' -v ss_col="$sample_num" -v chr_col="$chr_num" '
            NR==1 {next}
            $chr_col == "19" && $ss_col != "" && $ss_col ~ /^[0-9.]+$/ {vals[++c]=$ss_col}
            END {
                if (c > 0) {
                    n = asort(vals)
                    if (n % 2) {
                        med=vals[int(n/2)+1]
                    } else {
                        med=(vals[int(n/2)]+vals[int(n/2)+1])/2
                    }
                    low = med * 0.9
                    high = med * 1.1
                    print low, high
                }
            }
        '
    done \
    | awk '
        BEGIN {highest_min=0; lowest_max=""}
        {
            low=$1
            high=$2
            if(low > highest_min) highest_min = low
            if(lowest_max == "" || high < lowest_max) lowest_max = high
        }
        END {print highest_min, lowest_max}
    '
)

echo "highest_min: $highest_min"
echo "lowest_max: $lowest_max"

# Step 2: Build SNP table with extra flag column
output="chr19_snp_table.tsv"

{
    echo -e "snp\tntraits\tminmaf\tminss\tmaxss\tss_range_flag\ttrait_flag\tmaf_flag"

    awk -F',' '
    NR==1 {
        for(i=1; i<=NF; i++) {
            if($i=="raw_data_path") file_col=i;
            if($i=="snp") snp_col=i;
            if($i=="af") maf_col=i;
            if($i=="chrom") chrom_col=i;
            if($i=="sample_size") ss_col=i;
        }
        next;
    }
    NR>1 {
        print $file_col, $snp_col, $maf_col, $chrom_col, $ss_col;
    }
    ' "$gwas_info_file" \
    | while read trait_file snp_col_name maf_col_name chrom_col_name ss_col_name; do

        [[ "$trait_file" == *.gz ]] && cat_cmd="zcat" || cat_cmd="cat"

        header=$($cat_cmd "$trait_file" | head -1)
        IFS=$'\t' read -r -a cols <<< "$header"

        snp_col_idx=""
        maf_col_idx=""
        chrom_col_idx=""
        ss_col_idx=""
        for i in "${!cols[@]}"; do
            [[ "${cols[$i]}" == "$snp_col_name" ]] && snp_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$maf_col_name" ]] && maf_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$chrom_col_name" ]] && chrom_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$ss_col_name" ]] && ss_col_idx=$((i+1))
        done

        trait=$(basename "$trait_file" | sed 's/\..*//')

        $cat_cmd "$trait_file" | awk -F'\t' -v trait="$trait" \
            -v snp_col="$snp_col_idx" -v maf_col="$maf_col_idx" \
            -v chrom_col="$chrom_col_idx" -v ss_col="$ss_col_idx" \
            'NR > 1 && $chrom_col == "19" { print $snp_col "\t"trait"\t"$maf_col"\t"$ss_col }'

    done \
    | awk -F'\t' -v highest_min="$highest_min" -v lowest_max="$lowest_max" '
    {
      snp=$1
      trait=$2
      maf=$3
      ss=$4
      traits[snp][trait]=1
      if( (snp in minmaf) == 0 || maf < minmaf[snp]) minmaf[snp] = maf
      if( (snp in minss) == 0 || ss < minss[snp]) minss[snp] = ss
      if( (snp in maxss) == 0 || ss > maxss[snp]) maxss[snp] = ss
    }
    END {
      for(snp in traits) {
        ntraits = 0
        for(trait in traits[snp]) ntraits++
        flag = ((minss[snp] >= highest_min && maxss[snp] <= lowest_max) ? 1 : 0)
        traitflag = (ntraits == 3 ? 1 : 0)
        maf_flag = (minmaf[snp] > 0.05 ? 1 : 0)
        printf "%s\t%s\t%.5f\t%s\t%s\t%s\t%s\t%s\n", snp, ntraits, minmaf[snp], minss[snp], maxss[snp], flag, traitflag, maf_flag
      }
    }
    '
} > "$output"
