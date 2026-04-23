library(dplyr)
library(ieugwasr)

library(readr)

#X <- readRDS(snakemake@input[["zmat"]])
#X <- read.table(snakemake@input[["zmat"]], header = TRUE, sep = "\t")
#r2_thresh <- as.numeric(snakemake@wildcards[["r2_thresh"]])
#clump_kb <- snakemake@wildcards[["kb"]]
#ref_path  <- snakemake@params[["ref_path"]]
#out <- snakemake@output[["out"]]
#p <- snakemake@wildcards[["p"]]
#pthresh <- as.numeric(snakemake@params[["pthresh"]])
#is_mvmr <- as.numeric(snakemake@params[["is_mvmr"]])

# --- temp workdir for testing cleanliness --
workdir <- paste0("2_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

snps_in_ref <- file.path(workdir, "chr1_8traits_snps_in_ref_bim_file.tsv")
ld_pruned_out <- file.path(workdir, "ld_pruned_results.RDS")

# --- get just SNPs present in all traits AND ld ref files ---
print('filtering universal_snps.txt to only those found in ld ref files')

# gather inputs
# I had to add from ld pruning options:
ref_path <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal.bim"

# All lines from ld_file.tsv where the second column matches a value in the first column of snps_pass_all_filts.txt.

# other options
r2_thresh = 0.01 # list ok
clump_kb = 1000 # list ok
p = "pvalue" # list ok,  "pvalue" or "rank"
pthresh = 1
#ref_path: "1kg_plink/EUR" # list not ok
is_mvmr = 0

# --- helpful funcs, should be in another file 
# taken from GFA

# ---
if(!p %in% c("pvalue", "rank")){
  stop("Unknown prioritization option.\n")
}

gwas_info <- read_csv("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
traits <- gwas_info$name

# eventually change so pass_all_filt has 2 cols
tsv <- "../bash/keep_1_workdir_20260417_122646/final_pass_snps.tsv"
bim <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal.bim"
out <- "snps_pass_all_filts_in_ref.txt"

awk_prog <- 'BEGIN{OFS="\\t"}
NR==FNR { key=$1; row[key]=$0; next }
($2 in row) { keep[$2]=1 }
END { for (k in keep) print row[k] }'

get_snps_and_Z <- sprintf("awk '%s' %s %s > %s",
               awk_prog,
               shQuote(tsv), shQuote(bim), shQuote(out))

system(get_snps_and_Z)


#if(is_mvmr == 1){
#  Z_hat <- Z_hat[,-1,drop = FALSE]
#}

snps_and_Z <- read.table('snps_pass_all_filts_in_ref.txt')

if(p == "pvalue"){
  myp <- 2*pnorm(-abs(snps_and_Z[,2]))

#}else if(p == "rank"){
#  Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
#  min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
#  myp <- min_rank/max(min_rank)

}

dat <- data.frame(rsid = snps_and_Z[,1], pval = myp)

ref_pref <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal"
dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_pref)

pruned_snps <- snps_and_Z[which(snps_and_Z[,1] %in% dat_clump$rsid),1]

print(pruned_snps)

#saveRDS(X, file=out)

