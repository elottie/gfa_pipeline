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
workdir <- paste0("2_ld_prune_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

snps_in_ref <- file.path(workdir, "chr1_8traits_snps_in_ref_bim_file.tsv")
ld_pruned_out <- file.path(workdir, "ld_pruned_results.RDS")

# --- get just SNPs present in all traits AND ld ref files ---
print('filtering universal_snps.txt to only those found in ld ref files')

# gather inputs
# I had to add from ld pruning options:
ref_path <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal.bim"

# All lines from ld_file.tsv where the second column matches a value in the first column of snps_pass_all_filts.txt.
# I think i could take this away and put it in process_trait
awk_snps_in_ref <- paste(
  "cat", ref_path, "|",
  "awk -F'\\t' 'NR==FNR {snps[$1]=1; next} ( $2 in snps ) {print $2}'",
  "../bash/1_comb_form_chr1_8traits/snps_pass_all_filts.txt -",
  "| sort -u",
  ">", snps_in_ref
)

system(awk_snps_in_ref)

# other options
r2_thresh = 0.01 # list ok
clump_kb = 1000 # list ok
p = "pvalue" # list ok,  "pvalue" or "rank"
pthresh = 1
#ref_path: "1kg_plink/EUR" # list not ok
is_mvmr = 0

# --- helpful funcs, should be in another file 
# taken from GFA
align_beta <- function(X, beta_hat_name, af_name, A1_name, A2_name, upper=TRUE){
  A1 = X[[A1_name]]
  A2 = X[[A2_name]]
	
  flp = c("A" = "T", "G" = "C", "T" = "A",
          "C" = "G", "a"  = "t", "t" = "a",
          "c" = "g", "g" = "c")
  if(upper){
    X <- X %>% mutate( flip_strand = A1 == "T" | A2 == "T")
  }else{
    X <- X %>% mutate( flip_strand = A1 == "t" | A2 == "t")
  }

  if(missing(af_name)){
    X <- mutate(X, af = NA)
    af_name <- "tempaf"
    af_missing <- TRUE
  }else{
    af_missing <- FALSE
  }


  #print(head(X))
  X <- X %>% mutate(A1flp = case_when(flip_strand ~ flp[A1],
                                      TRUE ~ A1),
                    A2flp = case_when(flip_strand ~ flp[A2],
                                      TRUE ~ A2),
                    # afflp = case_when(flip_strand ~ 1-get(af_name),
                    #                    TRUE ~ get(af_name)),
                    tempbh = case_when(A1flp == "A" | A1flp == "a" ~ get(beta_hat_name),
                                     TRUE ~ -1*get(beta_hat_name)),
                    tempaf = case_when(A1flp == "A" | A1flp == "a" ~ get(af_name),
                                       TRUE ~ 1-get(af_name))) %>%
    select(-all_of(c(A1_name, A2_name))) %>%
    select(-all_of(c(af_name, beta_hat_name))) %>%
    mutate(A1 = case_when(A1flp == "A" | A1flp=="a" ~ A1flp,
                          TRUE ~ A2flp),
           A2 = case_when(A1flp == "A" | A1flp=="a" ~ A2flp,
                          TRUE ~ A1flp)) %>%
   select(-A1flp, -A2flp, -flip_strand)
  print(head(X)) 

  ix <- which(names(X)== "tempbh")
  names(X)[ix] <- beta_hat_name
  ix <- which(names(X)== "tempaf")
  names(X)[ix] <- af_name
  if(af_missing) X <- select(X, -tempaf)
  return(X)
}
process_trait <- function(trait, snps_in_ref, 
                          base_dir = "../../../gwas_summary_statistics/METSIM/with_rsid") {

  full_trait <- file.path(base_dir, trait, paste0(trait, "_regenie_rsid.tsv.gz"))
  message("... processing trait: ", full_trait)

  cmd <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || snps[$1]' %s -",
    shQuote(full_trait), shQuote(snps_in_ref)
  )

  filtered_trait <- read.table(pipe(cmd), header = TRUE, sep = "\t")
  print('head filtered_trait fter reading in:')
  print(head(filtered_trait))

  filtered_trait_harmon <- align_beta(filtered_trait,"BETA","A1FREQ","ALLELE1","ALLELE0")
  print('head filtered_trait_harmon:')
  print(head(filtered_trait_harmon))

  filtered_trait <- filtered_trait_harmon

  filtered_trait$Z <- filtered_trait$BETA / filtered_trait$SE
  print('head filtered_trait with z:')
  print(head(filtered_trait))

  list(
    A1 = filtered_trait$ALLELE1,
    A2 = filtered_trait$ALLELE0,
    Z = filtered_trait$Z,
    se = filtered_trait$SE,
    af = filtered_trait$A1FREQ,
    med_ss = median(filtered_trait$N)
  )

}


# ---
if(!p %in% c("pvalue", "rank")){
  stop("Unknown prioritization option.\n")
}

gwas_info <- read_csv("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
traits <- gwas_info$name

results <- list()
#for (trait in traits){
for (i in 1:1){
    trait = traits[i]
    results[[trait]] <- process_trait(trait, snps_in_ref)
}

#Z_hat <- X %>%
#  select(ends_with(".z")) %>%
#  as.matrix()
Z_hat <- as.data.frame(lapply(results, function(x) x$Z))

#if(is_mvmr == 1){
#  Z_hat <- Z_hat[,-1,drop = FALSE]
#}

if(p == "pvalue"){
  zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
  myp <- 2*pnorm(-abs(zmax))

#}else if(p == "rank"){
#  Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
#  min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
#  myp <- min_rank/max(min_rank)

}

snps <- scan(snps_in_ref, what = "", quiet = TRUE)
#dat <- data.frame(rsid = X$snp, pval = myp)
dat <- data.frame(rsid = snps, pval = myp)

ref_pref <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal"
dat_clump <- ld_clump(dat = dat,
                     clump_r2 = r2_thresh,
                     clump_p = pthresh,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_pref)

#ix <- which(X$snp %in% dat_clump$rsid)
#X <- X[ix,]

ix <- which(snps %in% dat_clump$rsid)

A1 <- results[[1]]$A1
A2 <- results[[1]]$A2

test <- cbind(
  data.frame(snp = snps[ix], A1 = A1[ix], A2 = A2[ix]),
  Z_hat[ix, , drop = FALSE]
)   # keeps Z columns numeric

print(head(test))
print(length(Z_hat))
print(length(test))

#saveRDS(X, file=out)

