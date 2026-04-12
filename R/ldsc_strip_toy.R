library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(GFA)
# library(stringr)

# added
library(data.table)
library(bench)
library(jsonlite)

# --- for once I graduate to snakemake ---
# sample size affects genetic covariance and h2 but not intercept or genetic correlation
gwas_info <- read_csv(snakemake@input[["gwas_info"]])
strip_list <- readRDS(snakemake@input[["strip_file"]])
strip_number <- as.numeric(snakemake@wildcard[["strip_num"]])
out <- snakemake@output[["out"]]
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])

# --- snakemake input ---
#rule R_ldsc_full:
#    input: Z = expand(data_dir + "{{prefix}}_zmat.{chrom}.RDS", chrom = range(1, 23)),
#           gwas_info = info_input,
#           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
#           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    params: cond_num = cond_num
#    wildcard_constraints: pt = r"[\d.]+"
#    script: "R/3_R_ldsc_all.R"

# --- function for processing traits ---
map_gwas_info_cols <- function(gwas_info, trait, trait_col = "name") {
  row <- gwas_info[gwas_info[[trait_col]] == trait, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected 1 row for trait=", trait, " found ", nrow(row))

  as.list(row[1, c("beta_hat","se","A1","A2","af")])
}

# Flip signs and strands so that allele 1 is always A (modifies X in place)
align_beta_dt <- function(X, upper = TRUE,
                          beta_col = "beta_hat",
                          af_col   = "af") {
  stopifnot(data.table::is.data.table(X))
  stopifnot(all(c("A1", "A2") %in% names(X)))
  stopifnot(beta_col %in% names(X))

  flp <- c("A"="T","G"="C","T"="A","C"="G",
           "a"="t","t"="a","c"="g","g"="c")

  af_present <- af_col %in% names(X)
  if (!af_present) {
    X[, (af_col) := NA_real_]  # create temporarily so code can be uniform
  }

  X[, ..__alignbeta_flip__ :=
       if (upper) (A1 == "T" | A2 == "T") else (A1 == "t" | A2 == "t")]

  X[, `:=`(
    ..__alignbeta_A1flp__ = data.table::fifelse(..__alignbeta_flip__, flp[A1], A1),
    ..__alignbeta_A2flp__ = data.table::fifelse(..__alignbeta_flip__, flp[A2], A2)
  )]

  X[, ..__alignbeta_condA__ := ..__alignbeta_A1flp__ %chin% c("A","a")]

  # beta_hat
  X[, (beta_col) := {
    b <- .SD[[1L]]
    if (!is.numeric(b)) b <- as.numeric(b)
    data.table::fifelse(..__alignbeta_condA__, b, -b)
  }, .SDcols = beta_col]

  # af
  X[, (af_col) := {
    p <- .SD[[1L]]
    if (!is.numeric(p)) p <- as.numeric(p)
    data.table::fifelse(..__alignbeta_condA__, p, 1 - p)
  }, .SDcols = af_col]

  # swap alleles to make A1 the A allele (after flip)
  X[, `:=`(
    A1 = data.table::fifelse(..__alignbeta_condA__, ..__alignbeta_A1flp__, ..__alignbeta_A2flp__),
    A2 = data.table::fifelse(..__alignbeta_condA__, ..__alignbeta_A2flp__, ..__alignbeta_A1flp__)
  )]

  X[, c("..__alignbeta_flip__", "..__alignbeta_A1flp__", "..__alignbeta_A2flp__", "..__alignbeta_condA__") := NULL]

  if (!af_present) X[, (af_col) := NULL]  # remove if it wasn't originally there

  invisible(X)
}

get_harmon_z <- function(gwas_info, trait, snps_in_ref, base_dir = "../../../gwas_summary_statistics/METSIM/with_rsid") {

  full_trait <- file.path(base_dir, trait, paste0(trait, "_regenie_rsid.tsv.gz"))
  message("... processing trait: ", full_trait)

  cmd <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || snps[$1]' %s -",
    shQuote(full_trait), shQuote(snps_in_ref)
  )

  m <- map_gwas_info_cols(gwas_info, trait)
  filtered_trait <- fread(cmd = cmd, sep = "\t", header = TRUE,
            select = unname(unlist(m)))  # pick what you need

  data.table::setnames(filtered_trait, old = unname(unlist(m)), new = names(m))

  print('head filtered_trait after reading in:')
  print(head(filtered_trait))
 
#  filtered_trait_harmon <- GFA:::align_beta(filtered_trait,beta_name,af_name) 
  filtered_trait_harmon <- align_beta_dt(filtered_trait) 
  print('head filtered_trait_harmon:')
  print(head(filtered_trait_harmon))

  filtered_trait <- filtered_trait_harmon
 
  # ought to guard for beta_hat & se anomalies, set to NA to preserve order 
  filtered_trait[, Z := beta_hat / se]
  print('head filtered_trait with z:')
  print(head(filtered_trait))

  return(filtered_trait[["Z"]])
}

# --- temp workdir for testing cleanliness --
workdir <- paste0("3_ldsc_strip_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

snps_in_ref <- file.path(workdir, "chr1_8traits_snps_in_ld_file.tsv")
ldsc_out <- file.path(workdir, "ldsc_results.RDS")

# --- get just SNPs present in all traits AND ld ref files ---
print('filtering universal_snps.txt to only those found in ld ref files')

# gather inputs
l2_dir <- "/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/"
chroms <- 1
m_files <- paste0(l2_dir, chroms, ".l2.M_5_50")
ld_files <- paste0(l2_dir, chroms, ".l2.ldscore.gz")

# All lines from ld_file.tsv where the second column matches a value in the first column of snps_pass_all_filts.txt.
awk_snps_in_ref <- paste(
  "zcat", ld_files, "|",
  "awk -F'\\t' 'BEGIN{OFS=\"\\t\"; print \"snp\", \"l2\"} NR==FNR {snps[$1]=1; next} snps[$2] {print $2, $6}' ../bash/1_comb_form_chr1_8traits/snps_pass_all_filts.txt -",
  ">", snps_in_ref
)

system(awk_snps_in_ref)

#l2 is col 2

# if M is num of variants used to compute ld scores, should be constant?
M <- purrr:::map(1, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

print('completed ld filtering')

# --- set up blocks for input to ldsc_strip.  use gwas info file ---

# read in gwas_info for trait names
gwas_info <- fread("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
traits <- gwas_info$name

# read in sample size stats
# eventually just get column 1 and 2
trait_ss <- fread("../bash/1_comb_form_chr1_8traits/trait_sample_stats.tsv",select=c(1,3))
setnames(trait_ss, c("trait", "ss_median"))
setkey(trait_ss, trait)   # enables fast lookup

# pass in strip_list
strip_list <- fromJSON("ldsc_trait_sets.json", simplifyVector = FALSE)  # list of character vectors

print(paste("Received smart-selected nblocks:", length(strip_list)))
print('These blocks are:')
print(strip_list)

# one day snakemake wildcard
strip_number <- 1

# --- strip processing! ---
n_snps <- as.integer(system(paste0("wc -l < ",snps_in_ref), intern = TRUE)) - 1L   # number of SNPs after filtering, -1 for header
print(n_snps)
n_traits <- length(traits)

l2 <- as.numeric(scan(pipe(sprintf("awk -F'\t' 'NR>1{print $2}' %s", snps_in_ref)), what="character"))
print('l2:')
print(head(l2))

ldsc_results <- list()
#for(s2 in strip_number:nblocks){
#for(s2 in strip_number:2){
for(s2 in strip_number:1){
    #block1_traits <- NULL
    block2_traits <- NULL
	
    if(s2 == strip_number){
        # read data for set 1 (first block)
        block1_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 1 (block", s2, "), has traits:",
        paste(block1_traits, collapse=", "), "\n")

	Z_hat_b1 <- matrix(NA_real_, nrow = n_snps, ncol = n_traits,
                dimnames = list(NULL, traits))   # or list(snps, traits) if you want

        for (trait in block1_traits) {
	  z <- get_harmon_z(gwas_info, trait, snps_in_ref)  # must return length M in correct order
          stopifnot(length(z) == n_snps)
          Z_hat_b1[, trait] <- z  
        }
	
	# for ldsc:  add within-block comparisons for block1_traits
	comparisons <- as.data.frame(t(combn(block1_traits, 2)), stringsAsFactors = FALSE)
        names(comparisons) <- c("trait1", "trait2")
        cat("Unique within-block comparisons for block", s2, ":\n")
        print(comparisons)
        # Get unique trait codes and turn string trait names to numbers
        all_traits <- unique(as.character(block1_traits))
    }

    if(s2 > strip_number + 1){
        # read set 2 (second block)
        block2_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 2 (block", s2, "), has traits:",
        paste(block2_traits, collapse=", "), "\n")

        Z_hat_b2 <- matrix(NA_real_, nrow = n_snps, ncol = n_traits,
                dimnames = list(NULL, traits))   # or list(snps, traits) if you want


        for (trait in block2_traits) {
          z <- get_harmon_z(gwas_info, trait, snps_in_ref)  # must return length M in correct order
          stopifnot(length(z) == n_snps)
          Z_hat_b2[, trait] <- z             
        }

	# for ldsc, pairwise comparison if both blocks are defined:
        comparisons <- expand.grid(trait1 = block1_traits,
                                   trait2 = block2_traits,
                                   stringsAsFactors = FALSE)
        cat("Pairwise comparisons between block", strip_number, "and", s2, ":\n")
        print(comparisons)
        # Get unique trait codes and turn string trait names to numbers
        all_traits <- unique(c(as.character(block1_traits), as.character(block2_traits)))

    }
    
    
    # Get unique trait codes and turn string trait names to numbers
    trait_to_idx <- setNames(seq_along(all_traits), all_traits)
    num_compar <- data.frame(
        trait1 = match(as.character(comparisons$trait1), all_traits),
        trait2 = match(as.character(comparisons$trait2), all_traits),
        row.names = NULL
    )
    print('Numeric comparisons:')
    print(num_compar)
    # easier to use dataframe later
    trait_map <- data.frame(
        idx = seq_along(all_traits),
        trait = all_traits,
        stringsAsFactors = FALSE
    )

    # need N, median of samples sizes, and Zs
    relev_Z <- Z_hat[, all_traits, drop = FALSE]
    print(head(relev_Z))
    # fast keyed lookup; returns in the same order as trait_order
    med_ss_vec <- trait_ss[colnames(relev_Z), ss_median]
    print(head(med_ss_vec))

    ldsc_result <- R_ldsc(
        Z_hat = relev_Z,
        ldscores = l2,
        N = med_ss_vec,
        ld_size = M,
	make_well_conditioned = FALSE,
        comparisons = num_compar
    )
    print(ldsc_result)

    # Store each result in your ldsc_results list
    block_comp_name <- paste0("block_", strip_number, "_vs_block_", s2)
    ldsc_results[[block_comp_name]] <- list("ldsc"=ldsc_result,"trait_name_map"=trait_map)
}

print(head(Z_hat))
str(ldsc_results)

saveRDS(ldsc_results, file=ldsc_out)



