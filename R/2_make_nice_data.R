library(data.table)

#rule make_nice_data:
#    input: gwas_info = info_input, 
#           pruned_snp_list = data_dir + "{prefix}_pruned_snps.{ldstring}.{chrom}.RDS"
#    params: usage = "gfa"  # would be MR for those which want beta & se
#    output: out = data_dir + "{prefix}_zmat.ldpruned_{ldstring}.{chrom}.RDS"
#    script: "R/make_nice_data.R"

gwas_info <- fread("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
pruned_snp_list = sprintf("../gfa_data/5e5Sig_Herit_Mets8_ld_pruned_chr%d.tsv", 1:22)
usage = "gfa"
out = paste0("../gfa_data/5e5Sig_Herit_Mets8_nice_data_for_",usage,".RDS")

# --- source helpful funcs ---
helper_path <- "harmon_helpers.R"
# eventually need to switch to this
#helper_path <- "R/harmon_helpers.R"
source(helper_path)

# --- temp workdir for testing cleanliness --
workdir <- paste0("2_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# temp output file defs
pruned_all_list <- file.path(workdir, "pruned_snps_across_chr.tsv")

# --- get pruned snps across chromsomes ---

# would not need the file.exists part, just for testing
pruned_all <- unique(unlist(lapply(pruned_snp_list[file.exists(pruned_snp_list)], readLines), use.names = FALSE))
writeLines(pruned_all, pruned_all_list)

# --- process ---
# read in data
#pruned_snps <- readRDS(pruned_snp_list)

#print(head(pruned_snps))

# handle the max snp thing? no, gfa does

# get Z and ss.  to get Z, need the harmon helper
traits <- gwas_info$name

Z_work <- matrix(NA_real_, length(pruned_all), length(traits),
                 dimnames = list(NULL, traits))
# run_gfa will be edited to accept this vector of medians so I don't need to write out whole ss since that's all it wants
ss_work <- matrix(NA_real_, 1, length(traits),
                 dimnames = list(NULL, traits))

for (trait in traits) {
  
  harmon <- harmon_dat(gwas_info, trait, pruned_all_list, return_ss=TRUE)

  Z_work[, trait] <- harmon$Z
  ss_work[, trait] <- median(harmon$ss)
}

# I don't believe we need to recreate this big input matrix with traitlabel.z and traitlabel.ss.  we can keep as 2 objects and reorder to match each other, better for mem
# well they are actually already ordered to match each other
# we just need to write them out!

head(Z_work)
head(ss_work)

if (identical(colnames(Z_work),colnames(ss_work))){
  save(Z_work, ss_work, file=out)
}

# --- delete temp workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
