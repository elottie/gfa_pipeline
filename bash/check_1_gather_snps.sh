#!/bin/sh
#SBATCH --account=jvmorr1

#zcat /nfs/turbo/sph-jvmorr/gwas_summary_statistics/METSIM/with_rsid/C999912707/C999912707_regenie_rsid.tsv.gz | awk -F'\t' '$2==14 {c++} END{print c}'
#zcat /nfs/turbo/sph-jvmorr/gwas_summary_statistics/METSIM/metsim_2022_10k/all_in/C999912707/C999912707_C999912707.regenie.gz | awk '$1==15 {c++} END{print c}'

csv="../5e5Sig_Herit_Mets_8ForLDSCStrip.csv"
k=15

pairs=$(mktemp)  # will hold: rsid<TAB>filepath  (unique per file later)
perfile=$(mktemp)

# 1) per-file counts + emit (rsid, file) pairs for rows with col2==15
tail -n +2 "$csv" | cut -d',' -f2 | while read -r f; do
  if [[ $f == *.gz ]]; then
    zcat "$f"
  else
    cat "$f"
  fi |
  awk -v k="$k" -v f="$f" -F'\t' '
    $2==k { c++; print $1 "\t" f }
    END   { print c+0 "\t" f > "/dev/stderr" }
  ' 2>>"$perfile" >>"$pairs"
done

echo "Per-file count of rows with col2==$k (count<TAB>file):"
sort -nr "$perfile"

# 2) across files: count in how many files each rsid appears with col2==15
#    (dedupe within-file so repeats in same file count once)
sort -u "$pairs" | awk -F'\t' '{n[$1]++} END{for (r in n) print r "\t" n[r]}' | sort -k2,2nr

rm -f "$pairs" "$perfile"
