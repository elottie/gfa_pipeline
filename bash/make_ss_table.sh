#!/bin/bash

gwas_info_file="../C100001554_And_Friends_3Metabolites.csv"
#head "$gwas_info_file"

# make the sample size table, +- 10% (eventually argument)
awk -F',' 'NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="raw_data_path") f_col=i;
        if($i=="sample_size") s_col=i;
        if($i=="chrom") c_col=i;
    }
    next
}
{
    print $f_col "\t" $s_col "\t" $c_col
}' "$gwas_info_file" \
| while IFS=$'\t' read file sample_col chrom_col; do
#    echo "file='$file', sample_col='$sample_col', chrom_col='$chrom_col'"
    trait=$(basename "$file" | sed 's/\..*//')

#    echo "$trait"
    if [[ "$file" == *.gz ]]; then
        cat_cmd="zcat"
    else
        cat_cmd="cat"
    fi
#    echo "$cat_cmd"

    # Grab column numbers for sample size and chromosome from the header row of the trait file:
    sample_num=$($cat_cmd "$file" | awk -F'\t' -v col="$sample_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')
    chr_num=$($cat_cmd "$file" | awk -F'\t' -v col="$chrom_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')

#    echo "$sample_num $chr_num"

    # Now filter for chromosome 19 and compute median sample size:
    $cat_cmd "$file" \
    | gawk -F'\t' -v col="$sample_num" -v chr_col="$chr_num" -v trait="$trait" '
    NR==1 {next}
    $chr_col == "19" && $col != "" && $col ~ /^[0-9.]+$/ {vals[++c]=$col}
    END {
        n = asort(vals)
        if (n == 0) exit
        if (n % 2) {
            med=vals[int(n/2)+1]
        } else {
            med=(vals[int(n/2)]+vals[int(n/2)+1])/2
        }
        low=med*0.9
        high=med*1.1
        printf "%s\t%.3f\t%.3f\t%.3f\n", trait, med, low, high
    }'
done > chr19_ss_table.tsv
