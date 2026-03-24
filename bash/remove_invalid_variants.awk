# USAGE:
# awk -F"\t" -v compute_pval="TRUE" -v snp_name="rsid" -v beta_name="beta" -v se_name="se" -v A1_name="A1" -v A2_name="A2" -v pos_name="pos" -v pval_name="p_value" -v ss_name="sample_size" -v af_name="allele_freq" -f gwas_format.awk input.tsv > output.tsv

BEGIN {
    OFS = "\t"
    # Accept col names via -v
    need_snp = (snp_name == "" ? "snp" : snp_name)
    need_A1 = (A1_name == "" ? "A1" : A1_name)
    need_A2 = (A2_name == "" ? "A2" : A2_name)
    need_beta = (beta_name == "" ? "beta_hat" : beta_name)
    need_se = (se_name == "" ? "se" : se_name)
#    need_pos = (pos_name == "" ? "pos" : pos_name)
#    need_pval = (pval_name == "" ? "p_value" : pval_name)
    need_ss = (ss_name == "" ? "sample_size" : ss_name)
    need_af = (af_name == "" ? "allele_freq" : af_name)
}

# --- Parse header and map colnames to indices ---
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[ $i ] = i
    }

    # validate required columns exist
    if (!(need_snp  in col)) { print "ERROR: missing column " need_snp  > "/dev/stderr"; exit 2 }
    if (!(need_A1   in col)) { print "ERROR: missing column " need_A1   > "/dev/stderr"; exit 2 }
    if (!(need_A2   in col)) { print "ERROR: missing column " need_A2   > "/dev/stderr"; exit 2 }
    if (!(need_beta in col)) { print "ERROR: missing column " need_beta > "/dev/stderr"; exit 2 }
    if (!(need_se   in col)) { print "ERROR: missing column " need_se   > "/dev/stderr"; exit 2 }
    if (!(need_ss   in col)) { print "ERROR: missing column " need_ss   > "/dev/stderr"; exit 2 }
    if (!(need_af   in col)) { print "ERROR: missing column " need_af   > "/dev/stderr"; exit 2 }

    # print the *current* header names for just these columns, in your chosen order
    print $(col[need_snp]), $(col[need_A1]), $(col[need_A2]), $(col[need_beta]), $(col[need_se]), $(col[need_ss]), $(col[need_af])
    next
}

# --- Main per-row logic ---
{
    # Access allele columns dynamically:
    snp  = $(col[need_snp])
    a1   = toupper($(col[need_A1]))
    a2   = toupper($(col[need_A2]))
    beta = $(col[need_beta])
    se   = $(col[need_se])
    ss   = $(col[need_ss])
    af   = $(col[need_af])

    # Remove illegal alleles (skip non-ATCG)
    if (!(a1 ~ /^[ACGT]$/ && a2 ~ /^[ACGT]$/)) next

    # Remove ambiguous alleles (AT, TA, GC, CG)
    if ( (a1=="A" && a2=="T") ||
         (a1=="T" && a2=="A") ||
         (a1=="C" && a2=="G") ||
         (a1=="G" && a2=="C") ) next

    # Remove duplicate SNPs
    # handled outside for lowest rss

    print snp, a1, a2, beta, se, ss, af
}

