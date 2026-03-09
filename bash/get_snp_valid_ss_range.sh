get_snp_valid_ss_range() {
    local gwas_info_file="$1"   # path to CSV file
    local trait_table="trait_ss_range_table.tsv"   # output file for table
    local highest_min lowest_max

    # write header to your table file
    echo -e "trait\tlow_tol_sample_size\tmedian_sample_size\thigh_tol_sample_size" > "$trait_table"

    # collect relevant columns
    awk -F',' 'NR==1 {
        for(i=1;i<=NF;i++) {
            if($i=="raw_data_path") f_col=i;
            if($i=="name") t_col=i;
            if($i=="sample_size") s_col=i;
            if($i=="chrom") c_col=i;
        }
        next
    }
    {print $f_col "\t" $t_col "\t" $s_col "\t" $c_col}
    ' "$gwas_info_file" |
while IFS=$'\t' read file trait_col sample_col chrom_col; do
        if [[ "$file" == *.gz ]]; then
            cat_cmd="zcat"
        else
            cat_cmd="cat"
        fi

        # Get column numbers for sample size and chromosome
        echo "$file"

        trait="$trait_col"
        sample_num=$($cat_cmd "$file" | awk -F'\t' -v col="$sample_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')
        chr_num=$($cat_cmd "$file" | awk -F'\t' -v col="$chrom_col" 'NR==1{for(i=1;i<=NF;i++) if($i==col) print i; exit}')

        echo "$trait"
        echo "$sample_num"
        echo "$chr_num"

        # Compute median for chromosome 19 sample sizes
        ss_range=$($cat_cmd "$file" | awk -F'\t' -v ss_col="$sample_num" -v chr_col="$chr_num" '
            BEGIN {OFS="\t"}
            NR==1{next}
            $chr_col == "19" && $ss_col != "" && $ss_col ~ /^[0-9.]+$/ {
                vals[++count] = $ss_col
            }
            END {
                if (count > 0) {
                    n = asort(vals)
                    if (n % 2) {
                        med=vals[int(n/2)+1]
                    } else {
                        med=(vals[int(n/2)]+vals[int(n/2)+1])/2
                    }
                    low = med * 0.9
                    high = med * 1.1
                    print low, med, high
                }
            }
        ')
        echo -e "${trait}\t${ss_range}" >> "$trait_table"

        done
    echo "Summary table written to $trait_table"

    # Now process the table
    highest_min_and_lowest_max=$(awk -F'\t' '
        NR==1 {next}
        {
            low = $2;
            high = $4;
            if (low > highest_min || NR==2) highest_min = low;
            if (lowest_max == "" || high < lowest_max) lowest_max = high;
        }
        END {print highest_min, lowest_max}
    ' "$trait_table")

    read highest_min lowest_max <<< "$highest_min_and_lowest_max"

    # Export results so they’re available to caller as variables
    printf "%s %s\n" "$highest_min" "$lowest_max"
}
