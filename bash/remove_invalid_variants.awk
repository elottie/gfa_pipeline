# USAGE:
# awk -F"\t" -v compute_pval="TRUE" -v snp_name="rsid" -v beta_name="beta" -v se_name="se" -v A1_name="A1" -v A2_name="A2" -v pos_name="pos" -v pval_name="p_value" -v ss_name="sample_size" -v af_name="allele_freq" -f gwas_format.awk input.tsv > output.tsv

BEGIN {
    OFS = "\t"
    # Accept col names via -v
    need_snp = (snp_name == "" ? "snp" : snp_name)
#    need_beta = (beta_name == "" ? "beta_hat" : beta_name)
#    need_se = (se_name == "" ? "se" : se_name)
    need_A1 = (A1_name == "" ? "A1" : A1_name)
    need_A2 = (A2_name == "" ? "A2" : A2_name)
#    need_pos = (pos_name == "" ? "pos" : pos_name)
#    need_pval = (pval_name == "" ? "p_value" : pval_name)
#    need_ss = (ss_name == "" ? "sample_size" : ss_name)
#    need_af = (af_name == "" ? "allele_freq" : af_name)
}

# --- Parse header and map colnames to indices ---
NR == 1 {
#    for (i = 1; i <= NF; i++) {
#        col[ $i ] = i
#    }
    # Optionally rename header columns
 #   $col[need_snp]   = "snp"
 #   $col[need_beta]  = "beta_hat"
 #   $col[need_se]    = "se"
 #   $col[need_A1]    = "A1"
 #   $col[need_A2]    = "A2"
 #   $col[need_pos]   = "pos"
 #   $col[need_pval]  = "p_value"
 #   $col[need_ss]    = "sample_size"
 #   $col[need_af]    = "allele_freq"

    # Output header
 #   print
    print "snp"
    next
}

# --- Main per-row logic ---
{
    # Access allele columns dynamically:
    a1 = $(col[need_A1])
    a2 = $(col[need_A2])
    snp = $(col[need_snp])

    # Uppercase alleles
    a1 = toupper(a1)
    a2 = toupper(a2)

    # Remove illegal alleles (skip non-ATCG)
    if (!(a1 ~ /^[ACGT]$/ && a2 ~ /^[ACGT]$/)) next

    # Remove ambiguous alleles (AT, TA, GC, CG)
    if ( (a1=="A" && a2=="T") ||
         (a1=="T" && a2=="A") ||
         (a1=="C" && a2=="G") ||
         (a1=="G" && a2=="C") ) next

    # Remove duplicate SNPs
    #snp_id = $col["snp"]
    if (snp_seen[snp]++) next

#    print
    print snp
}

