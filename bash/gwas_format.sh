#!/bin/bash

# Format GWAS summary statistics for CAUSE
# Usage example:

#source /path/to/gwas_format.sh

# Usage: gwas_format <input> <output> <snp> <beta_hat> <se> <A1> <A2> [chrom] [pos] [p_value] [sample_size] [allele_freq] [compute_pval]

#gwas_format "my_gwas.tsv" "out.tsv" "rsid" "beta" "se" "A1" "A2" "" "" "p_value" "sample_size" "allele_freq" "TRUE"

# debugging - - -
set -x

# source - - -
source get_col.sh

# -- Inner function: compute p-value --
compute_pvalue() {
    local infile="$1"
    local outfile="$2"

    awk -F'\t' '
    BEGIN{OFS="\t"} NR==1 {for(i=1;i<=NF;i++) idx[$i]=i; print; next}
    {
        # If p_value is missing and beta/se are present
        if ($idx["p_value"] == "" || $idx["p_value"] == "NA") {
            b = $idx["beta_hat"]
            s = $idx["se"]

            # Check if se is blank, NA, or 0
            if (s == "" || s == "NA" || s == 0) {
                $idx["p_value"] = "NA"
            } else {
                z = b / s
                p = 2 * (1 - systat(z))
                $idx["p_value"] = p
            }
        }
        print
    }
    function systat(z){return 0.5*erfc(-z/sqrt(2));}
    function erfc(x){return 1 - erf(x)}
    function erf(x){
      a1=0.254829592; a2=-0.284496736; a3=1.421413741; a4=-1.453152027; a5=1.061405429; p=0.3275911;
      sign=1; if(x<0) sign=-1; x=abs(x); t=1.0/(1.0+p*x);
      y=1.0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-x*x)
      return sign*y
    }
    ' "$infile" > "$outfile"
}

# -- Inner function: flip alleles --
flip_alleles() {
    local infile="$1"
    local outfile="$2"
    local a1_col="$3"
    local a2_col="$4"
    local beta_col="$5"
    local af_col="$6"

    awk -F'\t' -v a1_col="$a1_col" -v a2_col="$a2_col" -v beta_col="$beta_col" -v af_col="$af_col" -v OFS="\t" '
    function flip_allele(x) {
      return (x == "A") ? "T" : (x == "T") ? "A" : (x == "G") ? "C" : (x == "C") ? "G" : x
    }
    NR==1 {print $0; next}
    {
      A1 = $a1_col; A2 = $a2_col
      if (A1 == "T" || A2 == "T") {
        A1f = flip_allele(A1); A2f = flip_allele(A2)
      } else {
        A1f = A1; A2f = A2
      }
      if (A1f == "A") {
        eA1 = A1f; eA2 = A2f; new_beta = $beta_col
        if (af_col) new_af = ($af_col != "") ? $af_col : ""
      } else {
        eA1 = A2f; eA2 = A1f; new_beta = -$beta_col
        if (af_col) new_af = ($af_col != "") ? 1-$af_col : ""
      }
      for (i=1; i<=NF; i++) {
        if (i == a1_col) printf "%s", eA1
        else if (i == a2_col) printf "%s", eA2
        else if (i == beta_col) printf "%s", new_beta
        else if (af_col && i == af_col) printf "%s", new_af
        else printf "%s", $i
        printf (i == NF ? "\n" : OFS)
      }
    }
    ' "$infile" > "$outfile"
}

# -- Main workflow function --
gwas_format() {
    local input="$1"; shift
    local output="$1"; shift
    local SNP_COL="$1"; shift
    local BETA_COL="$1"; shift
    local SE_COL="$1"; shift
    local A1_COL="$1"; shift
    local A2_COL="$1"; shift
    local CHROM_COL="$1"; shift
    local POS_COL="$1"; shift
    local PVALUE_COL="$1"; shift
    local SAMPLESIZE_COL="$1"; shift
    local AF_COL="$1"; shift
    local compute_pval="$1"; shift

    TMP=$(mktemp)

    parse_header "$input"
    local snp_col=$(get_col "$SNP_COL")
    local beta_col=$(get_col "$BETA_COL")
    local se_col=$(get_col "$SE_COL")
    local a1_col=$(get_col "$A1_COL")
    local a2_col=$(get_col "$A2_COL")
    local chrom_col=$(get_col "$CHROM_COL")
    local pos_col=$(get_col "$POS_COL")
    local p_col=$(get_col "$PVALUE_COL")
    local ss_col=$(get_col "$SAMPLESIZE_COL")
    local af_col=$(get_col "$AF_COL")

    csvtk -t headers "$input"
    # 1. Renaming, uppercasing
    awk -v snp_col="$snp_col" \
    -v beta_col="$beta_col" \
    -v se_col="$se_col" \
    -v a1_col="$a1_col" \
    -v a2_col="$a2_col" \
    -v chrom_col="$chrom_col" \
    -v pos_col="$pos_col" \
    -v p_col="$p_col" \
    -v ss_col="$ss_col" \
    -v af_col="$af_col" \
    'NR==1{
      $snp_col="snp";
      $beta_col="beta_hat";
      $se_col="se";
      $a1_col="A1";
      $a2_col="A2";
      $chrom_col="chrom";
      $pos_col="pos";
      $p_col="p_value";
      $ss_col="sample_size";
      $af_col="allele_freq";
      print; next
    }
    {
      $a1_col = toupper($a1_col);
      $a2_col = toupper($a2_col);
      print;
    }' OFS='\t' "$input" > "$TMP"

    # 2. Remove illegal alleles
    awk -F'\t' -v a1_col="$a1_col" -v a2_col="$a2_col" '
    NR==1 {print; next}
    ($a1_col=="A" || $a1_col=="C" || $a1_col=="G" || $a1_col=="T") &&
    ($a2_col=="A" || $a2_col=="C" || $a2_col=="G" || $a2_col=="T")
    ' "$TMP" > "$TMP.legal"

    # 3. Remove ambiguous alleles
    awk -F'\t' -v a1_col="$a1_col" -v a2_col="$a2_col" '
    NR==1 {print; next}
    !( (($a1_col=="A" && $a2_col=="T") || ($a1_col=="T" && $a2_col=="A")) ||
       (($a1_col=="G" && $a2_col=="C") || ($a1_col=="C" && $a2_col=="G")) )
    {print}
    ' "$TMP.legal" > "$TMP.unambig"

    # 4. Remove duplicate SNPs
    awk -F'\t' -v snp_col="$snp_col" 'NR==1 {print; next} !a[$snp_col]++' "$TMP.unambig" > "$TMP.nodup"

    # 5. Calculate p-value (if requested)
    if [[ $compute_pval == "TRUE" ]]; then
      compute_pvalue "$TMP.nodup" "$TMP.pval"
      cp "$TMP.pval" "$TMP.nodup.step"
    else
      cp "$TMP.nodup" "$TMP.nodup.step"
    fi

    # 6. Flip alleles and effect
    flip_alleles "$TMP.nodup.step" "$output" "$a1_col" "$a2_col" "$beta_col" "$af_col"

    #rm -f "$TMP" "$TMP.legal" "$TMP.unambig" "$TMP.nodup" "$TMP.pval" "$TMP.nodup.step"
    echo "Formatted GWAS file written to $output"
}

# Main guard (call only if directly executed)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  gwas_format "$@"
fi
