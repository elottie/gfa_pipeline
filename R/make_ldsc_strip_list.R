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

gwas_info <- read.csv("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
out <- '../gfa_data/5e5Sig_Herit_Mets8_ldsc_strip_list.RDS'

# --- find out largest file size ---
files <- unique(na.omit(gwas_info$raw_data_path))
files <- files[file.exists(files) & grepl("\\.gz$", files)]
stopifnot(length(files) > 0)

find_file_sizes <- "xargs -d '\n' gzip -l | awk 'NR==1{next} $4==\"(totals)\" || $4==\"TOTAL\"{next} $2 ~ /^[0-9]+$/ {m=($2>m?$2:m)} END{print m}'"

max_uncomp <- as.numeric(system(find_file_sizes, intern = TRUE, input = paste(files, collapse = "\n"))) / 1e9
max_uncomp

# --- make trait sets ---
source('ldsc_strip_list_helpers.R')
# eventually need to switch to:
#source('R/ldsc_strip_list_helpers.R')

#"../C100001554_And_Friends_3Metabolites.csv"
#"../5e5Sig_Herit_Mets_8ForLDSCStrip.csv"
# set trait_memory_gb=1 for getting the 4 sets I want
sets <- make_trait_sets(gwas_info=gwas_info,memory_limit_gb=5,trait_memory_gb=max_uncomp)
str(sets)

saveRDS(sets,out)
