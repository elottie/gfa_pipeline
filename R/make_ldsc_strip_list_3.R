library(data.table)

#mem_limit_gb = 4
#rule make_trait_sets:
#    input: gwas_info = info_input
#    output: out = data_dir + "{prefix}_ldsc_strip_list.RDS"
#    params: mem_limit = mem_limit_gb
#    script: "R/make_ldsc_strip_list_3.R" 

gwas_info <- fread(snakemake@input[["gwas_info"]])
mem_limit <- as.numeric(snakemake@params[["mem_limit"]])
out <- snakemake@output[["out"]]

#gwas_info <- read.csv('../First8_Mets_ForLDSCStrip.csv')
#mem_limit <- 4
#out <- '../gfa_data/First8_Mets_ldsc_strip_list.RDS'

# --- make trait sets ---
source('R/ldsc_strip_list_helpers_3.R')
# eventually need to switch to:
#source('R/ldsc_strip_list_helpers.R')

#"../C100001554_And_Friends_3Metabolites.csv"
#"../5e5Sig_Herit_Mets_8ForLDSCStrip.csv"
# set trait_memory_gb=1 for getting the 4 sets I want
#in_gb <- 2688 / 1024
#in_gb <- 2750 / 1024
#in_gb <- 2800 / 1024
#in_gb <- 2915 / 1024
#in_gb <- 6000 / 1024

sets <- make_trait_sets(gwas_info=gwas_info,memory_limit_gb=mem_limit)
str(sets)

#sets_short <- list(sets[[1]][1:6])
#str(sets_short)
#saveRDS(sets_short,out)

saveRDS(sets,out)
