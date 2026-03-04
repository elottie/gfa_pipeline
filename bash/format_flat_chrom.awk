BEGIN {
    OFS = "\t"
    # Accept names for columns from -v
    seen_header = 0
}

# 1. Find column indices from header
NR == 1 {
    for (i = 1; i <= NF; i++) {
        colname = $i
        colnum[colname] = i
    }

    snp_col   = (colnum[snp_name]    ? colnum[snp_name]    : 0)
    chrom_col = (colnum[chrom_name]  ? colnum[chrom_name]  : 0)
    af_col    = (colnum[af_name]     ? colnum[af_name]     : 0)
    beta_col  = (colnum["BETA"]      ? colnum["BETA"]      : colnum["beta_hat"])
    # Prints header: If effect_is_or, add "beta" column
    if (effect_is_or == "TRUE" || effect_is_or == "true") {
        print $0, "beta"
    } else {
        print $0
    }
    next
}

{
    # 2. Filter rows to desired chromosome (if chrom_col and chrom are set)
    if (chrom_col && chrom != "" && $chrom_col != chrom) { next }

    # 3. Deduplication by SNP
    if (snp_col && $snp_col in unique_snps) { next }
    unique_snps[$snp_col] = 1

    # 4. AF filtering (if AF column and threshold provided)
    if (af_col && af_thresh != "") {
        afval = $af_col
        if (afval == "" || afval == "NA") next
        if (afval <= af_thresh || afval >= (1 - af_thresh)) next
    }

    # 5. Odds ratio log-transform if requested
    if (effect_is_or == "TRUE" || effect_is_or == "true") {
        # add log(beta_col) column
        log_beta = ($beta_col != "" && $beta_col != "NA") ? log($beta_col) : "NA"
        print $0, log_beta
    } else {
        print $0
    }
}
