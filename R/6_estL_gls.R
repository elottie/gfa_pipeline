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

# --- read in snakemake inputs and do checks ---
# snp_file is already dedup, headered, and ready to use from step 1
snps <- fread(snp_file, header = TRUE, select = 1)[[1]]

stopifnot(identical(rownames(R),colnames(R)))
nms_r <- colnames(R)
stopifnot(all(nms_r == gfafit$names))

# --- get Z and ss.  to get Z, need the harmon helper ---
source("R/harmon_helpers.R")
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

# subset to traits that made it through gfa
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
#res <- bind_cols(data[,c("chrom", "snp", "REF", "ALT")], res)

# here add back chrom, ref, alt
# a bit clunky, but chrom, ref, and alt do not change by trait.  so just take the last iteration of them from last harmon object
snps_dt <- data.table(chrom = harmon$chrom, snp = snps, ref = harmon$ref, alt = harmon$alt)
res <- bind_cols(snps_dt,res)
#saveRDS(res, file = out)
fwrite(res,file=out, sep="\t")

# --- delete temp workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
