#library(data.table)

# Example:
# blocks <- make_trait_sets("traits.csv", name_col = "name")
# str(blocks)

#def make_trait_sets(
#    csv_path: str,
#    name_col: str = "name",
#    memory_limit_gb: int = 150,
#    trait_memory_gb: int = 5,
#    blocks_at_once: int = 2,
#    min_traits_per_block: int = 2,
#) -> dict[str, list[str]]:

#rule make_trait_sets:
#    input: gwas_info = info_input
#    output: out = data_dir + "{prefix}_trait_sets.json"
#    params: mem_limit = mem_limit_gb
#    script: "R/make_ldsc_strip_list.py"

#find ../../../gwas_summary_statistics/METSIM/with_rsid -type f -name "*.tsv.gz" -exec gzip -l {} +
#	compressed        uncompressed  ratio uncompressed_name
#          402148080          1184878021  66.1% ../../../gwas_summary_statistics/METSIM/with_rsid/C100009051/C100009051_regenie_rsid.tsv
#          380906402          1129289232  66.3% ../../../gwas_summary_statistics/METSIM/with_rsid/C100002015/C100002015_regenie_rsid.tsv
#          402144894          1185053340  66.1% ../../../gwas_summary_statistics/METSIM/with_rsid/C100002397/C100002397_regenie_rsid.tsv

gwas_info <- read.csv('../First8_Mets_ForLDSCStrip.csv')
out <- '../gfa_data/First8_Mets_ldsc_strip_list.RDS'

# --- make trait sets ---
source('ldsc_strip_list_helpers_2.R')
# eventually need to switch to:
#source('R/ldsc_strip_list_helpers.R')

#"../C100001554_And_Friends_3Metabolites.csv"
#"../5e5Sig_Herit_Mets_8ForLDSCStrip.csv"
# set trait_memory_gb=1 for getting the 4 sets I want
#in_gb <- 2688 / 1024
#in_gb <- 2750 / 1024
in_gb <- 2800 / 1024
#in_gb <- 2915 / 1024
#in_gb <- 6000 / 1024
sets <- make_trait_sets(gwas_info=gwas_info,memory_limit_gb=in_gb)
str(sets)

#sets_short <- list(sets[[1]][1:6])
#str(sets_short)
#saveRDS(sets_short,out)

saveRDS(sets,out)
