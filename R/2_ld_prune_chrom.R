#library(dplyr)
library(ieugwasr)
#library(readr)

#X <- readRDS(snakemake@input[["zmat"]])
#X <- read.table(snakemake@input[["zmat"]], header = TRUE, sep = "\t")
#r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
#clump_kb <- snakemake@wildcards[["kb"]]
#ref_path  <- snakemake@params[["ref_path"]]
#out <- snakemake@output[["out"]]
#p <- snakemake@wildcards[["p"]]
#pthresh <- as.numeric(snakemake@params[["pthresh"]])
#is_mvmr <- as.numeric(snakemake@params[["is_mvmr"]])

snp_files <- sprintf("../gfa_data/First8_Mets_snps_chr%d.tsv", 1:22)
ref_path = "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal"
r2_thresh = 0.01 # list ok
clump_kb = 1000 # list ok
p = "pvalue" # list ok,  "pvalue" or "rank"
pthresh = 1
# is_mvmr goes away and is handled in step1
# this also goes away with snakemake --
chrom <- 1
# --
out <- paste0("../gfa_data/First8_Mets_ld_pruned_chr",chrom,".tsv")

# --- temp workdir for testing cleanliness --
workdir <- paste0("2_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

snps_in_bim <- paste0(workdir,"/snps_pass_all_filts_in_bim.txt")

# --- get snps that are in bim ---
bim <- paste0(ref_path,".bim")

awk_snps_in_bim <- 'BEGIN{OFS="\\t"}
  NR==FNR { key=$1; row[key]=$0; next }
  ($2 in row) { keep[$2]=1 }
  END { for (k in keep) print row[k] }'

get_snps_and_Z <- sprintf("cat %s | awk '%s' - %s > %s",
			  paste(shQuote(snp_files), collapse = " "),
			  awk_snps_in_bim, shQuote(bim), shQuote(snps_in_bim))
system(get_snps_and_Z)

snps_and_Z <- read.table(snps_in_bim)

# --- p-val stuff and ld clumping ---
if(p != "pvalue")){
  stop("Unaccepted prioritization option; only pvalue accepted.\n")
}

if(p == "pvalue"){
  myp <- 2*pnorm(-abs(snps_and_Z[,2]))

dat <- data.frame(rsid = snps_and_Z[,1], pval = myp)

dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)

pruned_snps <- snps_and_Z[which(snps_and_Z[,1] %in% dat_clump$rsid),1]

#print(head(pruned_snps))
#print(length(pruned_snps))

#saveRDS(pruned_snps, file=out)
writeLines(pruned_snps, out)

# --- delete temp workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
