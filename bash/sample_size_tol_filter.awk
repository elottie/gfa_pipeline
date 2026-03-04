# sample_size_tol_filter.awk
# Usage:
# awk -v ss_col=8 -v lower=5000 -v upper=8000 -f sample_size_tol_filter.awk infile.tsv > filtered.tsv

BEGIN {
    OFS = "\t"
}

NR==1 {
    print
    next
}

{
    ss = $ss_col
    # Only keep rows where sample_size is not NA and within bounds
    if (ss != "" && ss != "NA" && ss >= lower && ss <= upper)
        print
}
