library(dplyr)
library(purrr)
library(stringr)
library(GFA)


out <- snakemake@output[["out"]]
snp_file <- snakemake@input[["snp_list"]]
gfafit <- readRDS(snakemake@input[["gfa_file"]])
R <- readRDS(snakemake@input[["R"]])
stopifnot(all(R$names == gfafit$names))
R <- R$R

source("harmon_helpers.R")
# get Z and ss.  to get Z, need the harmon helper
traits <- gwas_info$name

snps <- readRDS(snp_file)

Z_hat <- matrix(NA_real_, length(snps), length(traits),
                 dimnames = list(snps, traits))

for (trait in traits) {
  harmon <- harmon_dat(gwas_info, trait, snp_file, return_ss=TRUE)
  Z_hat[, trait] <- harmon$Z
}

dat_names <- colnames(Z_hat)
ix <- which(dat_names %in% gfafit$names)
Z_hat <- Z_hat[,ix]
dat_names <- dat_names[ix]

o <- match(dat_names, gfafit$names)
Z_hat <- Z_hat[,o]
S <- matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat))
gls_sol <- GFA:::loadings_gls(Z_hat, S, R, gfafit$F_hat)
gls_zscores <- gls_sol$L/gls_sol$S
gls_pvals <- 2*pnorm(-abs(gls_zscores))
nf <- ncol(gfafit$F_hat)
res <- data.frame(cbind(gls_zscores, gls_pvals))
names(res) <- c(paste0("factor", 1:nf, ".z"), paste0("factor", 1:nf, ".p"))
res <- bind_cols(data[,c("chrom", "snp", "REF", "ALT")], res)
saveRDS(res, file = out)
