#!/bin/bash

# purpose:
# format GWAS summary statistics across all traits

# usage:
# source /path/to/gwas_format.sh
# gwas_format <input> <output> <snp> <beta_hat> <se> <A1> <A2> [chrom] [pos] [p_value] [sample_size] [allele_freq] [compute_pval]
# gwas_format "my_gwas.tsv" "out.tsv" "rsid" "beta" "se" "A1" "A2" "" "" "p_value" "sample_size" "allele_freq" "TRUE"

# --- debugging ---
set -x

# --- source helpers ---
source get_col.sh

# --- determine delimiter of gwas file input ---
gwas="$1"
gwas_delim=$(get_file_delimiter "$gwas")

# --- compute_p-value() --
# purpose:  compute 2-sided p-value for each SNP
# usage:  cat "$gwas" | compute_pvalue "$gwas_delim"
compute_pvalue() {
    local gwas_delim="$1"

    awk -v DELIM="$gwas_delim" '
      BEGIN {FS=DELIM; OFS=DELIM}
      # parse header and map column name to number
      NR==1 {for(i=1;i<=NF;i++) idx[$i]=i; print; next}
      # parse data
      {
        # if p_value is missing and beta/se are present, calculate p
        if ($idx["p_value"] == "" || $idx["p_value"] == "NA") {
            beta = $idx["beta_hat"]
            se = $idx["se"]
            # check if se is blank, NA, or 0
            if (se == "" || se == "NA" || se == 0) {
                $idx["p_value"] = "NA"
            } else {
                z = beta / se
                p = p_from_z(z)
                $idx["p_value"] = p
            }
        }
        print
    }
    # 2-sided p-val from z-score
    function p_from_z(z){return 0.5*norm_erf_c(-z/sqrt(2));}
    # complement of Abramowitz and Stegun Gaussian error function, for finding tail probabilities
    function norm_erf_c(x){return 1 - norm_erf(x)}
    # Abramowitz and Stegun Gaussian error function, normal_cdf = 1/2*[1+erf⁡(z/sqrt(2))]
    function norm_erf(x){
      # coefficients for approximation
      a1=0.254829592; a2=-0.284496736; a3=1.421413741; a4=-1.453152027; a5=1.061405429; p=0.3275911;
      sign = (x < 0) ? -1 : 1
      x=abs(x)
      t=1.0/(1.0 + p*x)
      y=1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-x*x)
    return sign*y
    }
    '
}

# -- flip_alleles() --
# purpose:  standardize ref and alt alleles and flip beta and af where necessary
# usage:  cat "$gwas" | flip_alleles "$gwas_delim" "$a1_col" "$a2_col" "$beta_col" "$af_col"
flip_alleles() {
    local gwas_delim="$1"
    local a1_col="$2"
    local a2_col="$3"
    local beta_col="$4"
    local af_col="$5"

    awk -v DELIM="$gwas_delim" -v a1="$a1_col" -v a2="$a2_col" -v beta="$beta_col" -v af="$af_col" '
      BEGIN {FS=DELIM; OFS=DELIM}
      function compl_allele(x) {
        return (x == "A") ? "T" : (x == "T") ? "A" : (x == "G") ? "C" : (x == "C") ? "G" : x
      }
      NR==1 {print $0; next}
      {
        A1 = $a1; A2 = $a2
        # flip strands
        if (A1 == "T" || A2 == "T") {
          A1_flip = compl_allele(A1); A2_flip = compl_allele(A2)
        } else {
          A1_flip = A1; A2_flip = A2
        }
        # we want all A1s to be A if possible
        if (A1_flip == "A") {
          new_A1 = A1_flip; new_A2 = A2_flip; new_beta = $beta
          # if the allele freq is something real (set):  condition ? value_if_true : value_if_false
          if (af) new_af = ($af != "") ? $af : "NA"
        # if A1 is not A after strand flip, flip alleles to try to get A1 to be A
        } else {
          new_A1 = A2_flip; new_A2 = A1_flip; new_beta = -$beta
          if (af) new_af = ($af != "") ? 1-$af : "NA"
        }
        # write out the new values for the columns where they need to go
        for (i=1; i<=NF; i++) {
          if (i == a1) printf "%s", new_A1
          else if (i == a2) printf "%s", new_A2
          else if (i == beta) printf "%s", new_beta
          else if (af && i == af) printf "%s", new_af
          else printf "%s", $i
          printf (i == NF ? "\n" : OFS)
        }
      }
   '
}

# -- gwas_format() --
# purpose:  main function to standardize gwas inputs, including ref and alt alleles
# usage:  gwas_format "${file_filt}_final" "$output" "$snp_name" "$beta_hat_name" "$se_name" "$A1_name" "$A2_name" "$chrom_name" "$pos_name" "$p_value_name" \
# "$sample_size_name" "$af_name" "TRUE"
gwas_format() {
    # read function inputs
    local gwas="$1"; shift
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

    # get column locations from names

if [[ ! -r "$file" ]]; then
        echo "Error: File '$file' does not exist or is not readable." >&2
        return 1
    fi

    if [[ "$file" == *.gz ]]; then
        header=$(awk 'NR==1 {print; exit}' <(gzip -cd "$file"))
    else
        header=$(awk 'NR==1 {print; exit}' "$file")
    fi


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
