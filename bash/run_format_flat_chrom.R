# Source the function definition
source("../R/format_ieu_chrom.R")

#/nfs/turbo/sph-jvmorr/gwas_summary_statistics/METSIM/with_rsid/C100001554/C100001554_regenie_rsid.tsv.gz \
#    1 0.05 \
#    RSID GENPOS CHROM ALLELE1 ALLELE0 BETA SE LOG10P A1FREQ N \
#    FALSE \
#    C100001554_format_flat_chrom.txt

# Call the function with your arguments
result <- format_flat_chrom(
  file = "/nfs/turbo/sph-jvmorr/gwas_summary_statistics/METSIM/with_rsid/C100001554/C100001554_regenie_rsid.tsv.gz",
  chrom = "1",
  af_thresh = 0.05,
  snp_name = "RSID",
  pos_name = "GENPOS",
  chrom_name = "CHROM",
  A1_name = "ALLELE1",
  A2_name = "ALLELE0",
  beta_hat_name = "BETA",
  se_name = "SE",
  p_value_name = "LOG10P",
  af_name = "A1FREQ",
  sample_size_name = "N",
  effect_is_or = FALSE
)

# Write the result to a file (CSV in this example)
write.csv(result, file = "C100001554_format_flat_chrom_Routput.csv", row.names = FALSE)

