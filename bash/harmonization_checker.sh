while read harmon_line; do
    harmon_key=$(echo "$harmon_line" | awk '{print $1}')
    R_output_line=$(awk -v val="$harmon_key" '$3 == val' C100001554_format_flat_chrom_Routput.csv)
    echo "$harmon_line" >> harmonization_check.txt
    echo "$R_output_line" >> harmonization_check.txt
done < harmonized_snps.txt
