# USAGE:
# awk -F"\t" -v compute_pval="TRUE" -v snp_name="rsid" -v beta_name="beta" -v se_name="se" -v A1_name="A1" -v A2_name="A2" -v pos_name="pos" -v pval_name="p_value" -v ss_name="sample_size" -v af_name="allele_freq" -f gwas_format.awk input.tsv > output.tsv

BEGIN {
    OFS = "\t"
    # Accept col names via -v
    need_snp = (snp_name == "" ? "snp" : snp_name)
    need_beta = (beta_name == "" ? "beta_hat" : beta_name)
    need_se = (se_name == "" ? "se" : se_name)
    need_A1 = (A1_name == "" ? "A1" : A1_name)
    need_A2 = (A2_name == "" ? "A2" : A2_name)
    need_pos = (pos_name == "" ? "pos" : pos_name)
    need_pval = (pval_name == "" ? "p_value" : pval_name)
    need_ss = (ss_name == "" ? "sample_size" : ss_name)
    need_af = (af_name == "" ? "allele_freq" : af_name)
}

# --- Parse header and map colnames to indices ---
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[ $i ] = i
    }
    # Optionally rename header columns
    $col[need_snp]   = "snp"
    $col[need_beta]  = "beta_hat"
    $col[need_se]    = "se"
    $col[need_A1]    = "A1"
    $col[need_A2]    = "A2"
    $col[need_pos]   = "pos"
    $col[need_pval]  = "p_value"
    $col[need_ss]    = "sample_size"
    $col[need_af]    = "allele_freq"

    # Output header
    print
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

    # Optionally compute p-value if requested
    if (compute_pval == "TRUE") {
        beta = $(col["need_beta"])
        se   = $(col["need_se"])
        pval = $(col["need_pval"])
        if (pval == "" || pval == "NA") {
            if (se != "" && se != "NA" && se != 0) {
                z = beta/se
                $col["need_pval"] = p_value_from_z(z)
            } else {
                $col["need_pval"] = "NA"
            }
        }
    }

    # Flip alleles and beta/AF if needed
    flip_alleles()

    print
}

# --------- Functions ---------
function compl_allele(x) {
    return (x == "A") ? "T" : (x == "T") ? "A" : (x == "G") ? "C" : (x == "C") ? "G" : x
}

function flip_alleles(   A1,A2,af,beta,A1_flip,A2_flip,new_A1,new_A2,new_beta,new_af) {
    A1 = $col["need_A1"]
    A2 = $col["need_A2"]
    beta = $col["need_beta"]
    af = $col["need_af"]

    # flip strands
    A1_flip = A1; A2_flip = A2
    if (A1 == "T" || A2 == "T") {
        A1_flip = compl_allele(A1); A2_flip = compl_allele(A2)
    }
    # we want all A1s to be A if possible
    if (A1_flip == "A") {
        new_A1 = A1_flip; new_A2 = A2_flip; new_beta = beta
        new_af = (af != "") ? af : "NA"
    } else {
        new_A1 = A2_flip; new_A2 = A1_flip; new_beta = -beta
        new_af = (af != "") ? 1-af : "NA"
    }A1 = $col["need_A1"]
    $col["need_A1"] = new_A1
    $col["need_A2"] = new_A2
    $col["need_beta"] = new_beta
    $col["need_af"] = new_af
}

function p_value_from_z(z) {
    return 2*(1 - normcdf(abs(z)))
}

function normcdf(x) { return 0.5 * (1 + erf(x / sqrt(2))) }
function erf(x){
    # Abramowitz and Stegun formula
    a1=0.254829592;a2=-0.284496736;a3=1.421413741
    a4=-1.453152027;a5=1.061405429;p=0.3275911
    sign = (x<0)?-1:1
    x = abs(x)
    t=1.0/(1.0+p*x)
    y=1.0 - ((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-x*x)
    return sign*y
}
