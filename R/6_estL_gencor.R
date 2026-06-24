library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)
library(data.table)

# also  needs to be edited bc of new z lists

# sample size affects genetic covariance and h2 but not intercept or genetic correlation
gwas_info_file <- snakemake@input[["gwas_info"]]
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

gwas_info <- fread(gwas_info_file)

# --- source helpful funcs ---
# make awks, sorts, and joins consistent across users
Sys.setenv(LC_ALL = "C")

harmon_helper_path <- "R/harmon_helpers.R"
ld_ref_helper_path <- "R/ld_ref_helpers.R"
# eventually need to switch to this
#helper_path <- "R/harmon_helpers.R"
source(harmon_helper_path)
source(ld_ref_helper_path)

# --- temp workdir for testing cleanliness -
workdir <- paste0("6_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# temp output file defs
snps_in_ref_file <- file.path(workdir, "snps_in_ld_file.tsv")

# --- get ld ref across chromsomes ---
concat_ld_ref <- build_concat_ld_ref(ld_files = ld_files,
                                     snps_in_ref_file = snps_in_ref_file,
                                     delim="\t")

concat_status <- system(concat_ld_ref)
if (concat_status != 0) {
  stop("Failed to create unsorted reference with header")
}

print('completed ld ref concatenation')

# --- other random input setup ---
# if M is num of variants used to compute ld scores, should be constant?
M <- purrr:::map(1, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

# read in gwas_info for trait names
traits <- gwas_info$name

# --- get Z and ss.  to get Z, need the harmon helper ---
# read in the ref snps for rownames
# I have already deduplicated them
snps_in_ref <- fread(snps_in_ref_file, header = TRUE, select = 1)[[1]]
n_snps <- length(snps_in_ref)
print(n_snps)
n_traits <- length(traits)

Z_hat <- matrix(NA_real_, n_snps, length(traits),
		                 dimnames = list(snps_in_ref, traits))
print(dim(Z_hat))

l2 <- as.numeric(scan(pipe(sprintf("awk 'NR > 1 {print $2}' %s", snps_in_ref_file)), what="character"))

for (trait in traits) {
  # don't need return ss, let it default to false
  # in make_nice_data and estL_gls we do not have need_invalid_snps_rm because the list we pass to harmon_dat is from snp lists:  already rmed invalid snps
  # here, as in the R step 1, we add need_invalid_snps_rm because the list we are passing in is snps that are in reference file
  # i.e. w/o removal optino, there are invalid snps that are found in trait (ex. G/C snp), found in ref file, which would be extracted and used in gencor, which we don't want
  harmon <- harmon_dat(gwas_info, trait, snps_in_ref_file, need_invalid_snps_rm=TRUE)

  # add check that snps are identical to rownames(Z_Hat)
  if (identical(harmon$snps,rownames(Z_hat))){
    Z_hat[, trait] <- harmon$Z
  } else {
    stop('rowname snps used in harmon_dat are not the same as rownames of destination Z_hat matrix')
  }
}


# --- obtain genetic correlation ---
R <- R_ldsc(Z_hat = Z_hat,
            ldscores = l2,
            ld_size = M,
            N = rep(1, ncol(Z_hat)),
            return_gencov = TRUE,
            make_well_conditioned = FALSE # not needed we only need Rg
)

# --- save output ---
saveRDS(R, file=out)

# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)

