library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)

# also  needs to be edited bc of new z lists

# sample size affects genetic covariance and h2 but not intercept or genetic correlation
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

ld <- purrr::map_dfr(1:22, function(c){
  read_table(ld_files[c])
})

M <- purrr:::map(1:22, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()


# --- get Z and ss.  to get Z, need the harmon helper ---
source("harmon_helpers.R")
#source("R/harmon_helpers.R")
traits <- gwas_info$name

Z_hat <- matrix(NA_real_, length(snps), length(traits),
		                 dimnames = list(snps, traits))

print(dim(Z_hat))
print(length(snps))
print(head(snps))

for (trait in traits) {

	  # don't need return ss, let it default to false
	  harmon <- harmon_dat(gwas_info, trait, snp_file, return_alleles=TRUE)

  # add check that snps are identical to rownames(Z_Hat)
  if (identical(harmon$snps,rownames(Z_hat))){
    Z_hat[, trait] <- harmon$Z
  } else {
    stop('rowname snps used in harmon_dat are not the same as rownames of destination Z_hat matrix')
  }
}

X <- map_dfr(z_files, function(f){
  readRDS(f) %>%
    rename(SNP = snp) %>%
    inner_join(., ld)})

Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nms <- str_replace(colnames(Z_hat), ".z$", "")

R <- R_ldsc(Z_hat = Z_hat,
            ldscores = X$L2,
            ld_size = M,
            N = rep(1, ncol(Z_hat)),
            return_gencov = TRUE,
            make_well_conditioned = FALSE # not needed we only need Rg
)

saveRDS(R, file=out)

