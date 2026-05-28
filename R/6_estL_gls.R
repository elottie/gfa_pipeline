library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(GFA)


out <- snakemake@output[["out"]]
gwas_info <- fread(snakemake@input[["gwas_info"]])
snp_file <- snakemake@input[["snp_list"]]
gfafit <- readRDS(snakemake@input[["gfa_file"]])
R <- readRDS(snakemake@input[["R"]])

# --- temp workdir for testing cleanliness --

workdir <- paste0("6_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

snp_file_no_head <- paste0(workdir,"/snps_no_head.txt")

# --- read in snakemake inputs and do checks ---
#snps <- fread(snp_file)[, 1, with=FALSE]
snps <- unique(fread(snp_file)[[1]])
writeLines(snps,snp_file_no_head)

stopifnot(identical(rownames(R),colnames(R)))
nms_r <- colnames(R)
stopifnot(all(nms_r == gfafit$names))

# --- get Z and ss.  to get Z, need the harmon helper ---
source("R/harmon_helpers.R")
traits <- gwas_info$name

Z_hat <- matrix(NA_real_, length(snps), length(traits),
                 dimnames = list(snps, traits))

print(dim(Z_hat))
print(length(snps))
print(head(snps))

for (trait in traits) {
  harmon <- harmon_dat(gwas_info, trait, snp_file_no_head, return_ss=TRUE)
  print(length(harmon$Z))
  print(head(harmon$Z))
  Z_hat[, trait] <- harmon$Z
}

nms_z <- colnames(Z_hat)
ix <- which(nms_z %in% gfafit$names)
Z_hat <- Z_hat[,ix]
nms_z <- nms_z[ix]

o <- match(nms_z, gfafit$names)  # in order of nms_z
Z_hat <- Z_hat[,o]

# --- gls stuff ---
S <- matrix(1, nrow = nrow(Z_hat), ncol = ncol(Z_hat))
gls_sol <- GFA:::loadings_gls(Z_hat, S, R, gfafit$F_hat)
gls_zscores <- gls_sol$L/gls_sol$S
gls_pvals <- 2*pnorm(-abs(gls_zscores))
nf <- ncol(gfafit$F_hat)
res <- data.frame(cbind(gls_zscores, gls_pvals))
names(res) <- c(paste0("factor", 1:nf, ".z"), paste0("factor", 1:nf, ".p"))

# --- add back to some snp info: chrom, snp, ref, alt, gls_z, gls_p ---
# old way gave you more info, but how necessary?  could lookup after if you really want
#res <- bind_cols(data[,c("chrom", "snp", "REF", "ALT")], res)

# here add back chrom, ref, alt
snps_dt <- data.table(snp = snps)
res <- bind_cols(snps_dt,res)
saveRDS(res, file = out)

# --- delete temp workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
