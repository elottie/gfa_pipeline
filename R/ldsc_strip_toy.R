library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(GFA)
# library(stringr)

# added
library(data.table)
library(bench)

# --- for once I graduate to snakemake ---
# sample size affects genetic covariance and h2 but not intercept or genetic correlation
# gwas_info <- read_csv(snakemake@input[["gwas_info"]])
# strip_list <- readRDS(snakemake@input[["strip_file"]])
# strip_number <- as.numeric(snakemake@wildcard[["strip_number"]])
# out <- snakemake@output[["out"]]
# ld_files <- unlist(snakemake@input[["l2"]])
# m_files <- unlist(snakemake@input[["m"]])

# --- make toy setup to test Jean's updated R_ldsc ---
# Toy Z-score matrix: rows = SNPs, columns = traits
#Z_hat <- matrix(c(
#  1.2, 0.8, 0.9, 1.9,  # SNP1
#  1.1, 0.7, 1.5, 0.2,  # SNP2
#  0.9, 0.5, 1.8, 0.7,  # SNP3
#  1.0, 0.6, 2.1, 0.9    # SNP4
#), nrow = 4, byrow = TRUE)
# LD scores (numeric vector) for each SNP.  come from LD ref panel, not me
#ldscores <- c(1.5, 2.0, 1.8, 2.2)
# LD size, num of variants used to compute LD scores, right now arbitrary scalar
#ld_size <- 4
#N: sample sizes per trait (if equal across SNPs, can be a vector)
#N <- c(1000, 1200, 1200, 1100)
#blocks <- NULL
#ncores <- 1
# the new option, comparisons of traits.  comparisons MUST be df, not matrix
#comparisons <- data.frame(trait1 = c(1,3),
#                          trait2 = c(2,4))
#print(comparisons)

#result <- R_ldsc(
#  Z_hat = Z_hat,
#  ldscores = ldscores,
#  ld_size = ld_size,
#  N = N,
#  blocks = blocks,
#  ncores = ncores,
#  make_well_conditioned = FALSE,
#  comparisons = comparisons
#)
#print(result)

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
get_col_map <- function(gwas_info, trait, trait_col = "name") {
  row <- gwas_info[gwas_info[[trait_col]] == trait, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected 1 row for trait=", trait, " found ", nrow(row))

  as.list(row[1, c("beta_hat","se","A1","A2","af")])
}

# align_beta <- function(X, a1_name, a2_name, beta_hat_name, af_name, upper=TRUE){
#   A1 = X[[a1_name]]
#   A2 = X[[a2_name]]
# 
#   flp = c("A" = "T", "G" = "C", "T" = "A",
#           "C" = "G", "a"  = "t", "t" = "a",
#           "c" = "g", "g" = "c")
#   if(upper){
#     X <- X %>% mutate( flip_strand = A1 == "T" | A2 == "T")
#   }else{
#     X <- X %>% mutate( flip_strand = A1 == "t" | A2 == "t")
#   }
# 
#   if(missing(af_name)){
#     X <- mutate(X, af = NA)
#     af_name <- "tempaf"
#     af_missing <- TRUE
#   }else{
#     af_missing <- FALSE
#   }
# 
# 
#   # CHANGED the select of A1 and A2
#   X <- X %>% mutate(A1flp = case_when(flip_strand ~ flp[A1],
#                                       TRUE ~ A1),
#                     A2flp = case_when(flip_strand ~ flp[A2],
#                                       TRUE ~ A2),
#                     # afflp = case_when(flip_strand ~ 1-get(af_name),
#                     #                    TRUE ~ get(af_name)),
#                     tempbh = case_when(A1flp == "A" | A1flp == "a" ~ get(beta_hat_name),
#                                      TRUE ~ -1*get(beta_hat_name)),
#                     tempaf = case_when(A1flp == "A" | A1flp == "a" ~ get(af_name),
#                                        TRUE ~ 1-get(af_name))) %>%
#     select(-a1_name, -a2_name) %>%
#     select(-all_of(c(af_name, beta_hat_name))) %>%
#     mutate(A1 = case_when(A1flp == "A" | A1flp=="a" ~ A1flp,
#                           TRUE ~ A2flp),
#            A2 = case_when(A1flp == "A" | A1flp=="a" ~ A2flp,
#                           TRUE ~ A1flp)) %>%
#     select(-A1flp, -A2flp, -flip_strand)
# 
#   ix <- which(names(X)== "tempbh")
#   names(X)[ix] <- beta_hat_name
#   ix <- which(names(X)== "tempaf")
#   names(X)[ix] <- af_name
#   if(af_missing) X <- select(X, -tempaf)
#   return(X)
# }
# 

# Flip signs and strands so that allele 1 is always A (modifies X in place)
align_beta_dt <- function(X, beta_hat_name, af_name, upper = TRUE) {
  stopifnot(is.data.table(X))
  stopifnot(all(c("A1","A2") %in% names(X)))
  stopifnot(beta_hat_name %in% names(X))

  flp <- c("A"="T","G"="C","T"="A","C"="G",
           "a"="t","t"="a","c"="g","g"="c")

  bh <- beta_hat_name
  af_missing <- missing(af_name)
  if (af_missing) {
    af <- "..__alignbeta_af__"
    X[, (af) := NA_real_]
  } else {
    stopifnot(af_name %in% names(X))
    af <- af_name
  }

  X[, ..__alignbeta_flip__ := if (upper) (A1 == "T" | A2 == "T") else (A1 == "t" | A2 == "t")]

  X[, `:=`(
    ..__alignbeta_A1flp__ = fifelse(..__alignbeta_flip__, flp[A1], A1),
    ..__alignbeta_A2flp__ = fifelse(..__alignbeta_flip__, flp[A2], A2)
  )]

  X[, ..__alignbeta_condA__ := ..__alignbeta_A1flp__ %chin% c("A","a")]

  # beta: be explicit and numeric
  X[, (bh) := {
    b <- .SD[[1L]]
    if (!is.numeric(b)) b <- as.numeric(b)
    fifelse(..__alignbeta_condA__, b, -b)
  }, .SDcols = bh]

  # af: similarly explicit
  X[, (af) := {
    p <- .SD[[1L]]
    if (!is.numeric(p)) p <- as.numeric(p)
    fifelse(..__alignbeta_condA__, p, 1 - p)
  }, .SDcols = af]

  X[, `:=`(
    A1 = fifelse(..__alignbeta_condA__, ..__alignbeta_A1flp__, ..__alignbeta_A2flp__),
    A2 = fifelse(..__alignbeta_condA__, ..__alignbeta_A2flp__, ..__alignbeta_A1flp__)
  )]

  X[, c("..__alignbeta_flip__", "..__alignbeta_A1flp__", "..__alignbeta_A2flp__", "..__alignbeta_condA__") := NULL]
  if (af_missing) X[, (af) := NULL]

  invisible(X)
}

process_trait <- function(trait, snps_in_ref,beta_name,se_name,af_name,gwas_info,
                          base_dir = "../../../gwas_summary_statistics/METSIM/with_rsid") {

  full_trait <- file.path(base_dir, trait, paste0(trait, "_regenie_rsid.tsv.gz"))
  message("... processing trait: ", full_trait)

  cmd <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || snps[$1]' %s -",
    shQuote(full_trait), shQuote(snps_in_ref)
  )

  m <- get_col_map(gwas_info, trait)
  filtered_trait <- fread(cmd = cmd, sep = "\t", header = TRUE,
            select = unname(unlist(m)))  # pick what you need

  data.table::setnames(filtered_trait, old = unname(unlist(m)), new = names(m))

  #filtered_trait <- read.table(pipe(cmd), header = TRUE, sep = "\t")
  print('head filtered_trait after reading in:')
  print(head(filtered_trait))
 
#  b <- bench::mark(
#  dplyr = {
#    X <- as_tibble(filtered_trait)   # or whatever your original input class is
#    GFA:::align_beta(X, beta_name, af_name)
#  },
#  dt = {
#    X <- data.table::copy(filtered_trait)
#    align_beta_dt(X, beta_name, af_name)
#    X
#  },
#  iterations = 15,
#  mem_alloc = TRUE,
#  check = FALSE
#  )

#  filtered_trait_harmon <- GFA:::align_beta(filtered_trait,beta_name,af_name) 
  filtered_trait_harmon <- align_beta_dt(filtered_trait,beta_name,af_name) 
  print('head filtered_trait_harmon:')
  print(head(filtered_trait_harmon))

  filtered_trait <- filtered_trait_harmon
  
  filtered_trait[, Z := beta_hat / se]
#  filtered_trait$Z <- filtered_trait[[beta_name]] / filtered_trait[[se_name]]
  print('head filtered_trait with z:')
  print(head(filtered_trait))

  list(
    Z = filtered_trait$Z#,
#    med_ss = NA #median(filtered_trait$N)
#     bench <- b
     )
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
# I had to add from ld pruning options:
bim_path <- "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink_hg38/EUR_autosomal.bim"

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
gwas_info <- read_csv("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
traits <- gwas_info$name

# set block size, one day based on our knowledge of memory taken by trait.  150 GB ish on std partition
# we want few blocks, but lower limit on blocks
# assuming each trait file is 5 GB unzipped (it will be smaller from awk wizardry) and we hold 2 blocks at a time, then how many traits can we have per block?
# 2 blocks * x traits / block * 5 gb / trait = 150 gb
# x = 15 traits / block maximum
# then we want to know how many blocks we can get away with
# b blocks min = n total traits / (15 traits / block)
# additional note:  each block needs at least 2 traits so we can make one pairwise comparison

memory_limit_gb <- 150     # total available memory
trait_memory_gb <- 5       # each trait file is 5 GB
blocks_at_once <- 2    # we hold two blocks at once

# Calculate max traits per block given memory constraint
max_traits_per_block <- floor(memory_limit_gb / (blocks_at_once * trait_memory_gb)) # = 15
min_traits_per_block <- 2

ntraits <- length(traits)

# Find minimum number of blocks so every block has at least min_traits_per_block AND not more than traits_per_block
min_nblocks <- ceiling(ntraits / max_traits_per_block)        # Minimum blocks, using largest block possible

# Pick the smallest nblocks possible according to memory constraint and minimum traits per block
#nblocks <- max(min_nblocks, 1) # avoid 0 blocks
# I INTERVENE bc I want four right now for testing:
nblocks <- 4

block_size <- ceiling(ntraits / nblocks)

# Assign traits to blocks
strip_indices <- rep(1:nblocks, each=block_size, length.out=ntraits)
strip_list <- split(traits, strip_indices)
names(strip_list) <- paste0("set", seq_along(strip_list))

# Guarantee no single-trait or empty block at the end
while(length(strip_list[[length(strip_list)]]) < min_traits_per_block && length(strip_list) > 1) {
  strip_list[[length(strip_list)-1]] <- c(strip_list[[length(strip_list)-1]], strip_list[[length(strip_list)]])
  strip_list <- strip_list[-length(strip_list)]
  names(strip_list) <- paste0("set", seq_along(strip_list))
}

print(paste("Smart-selected nblocks:", length(strip_list)))
print('These blocks are:')
print(strip_list)

# one day snakemake wildcard
strip_number <- 1

# --- strip processing! ---
results <- list()
ldsc_results <- list()
#for(s2 in strip_number:nblocks){
#for(s2 in strip_number:2){
for(s2 in strip_number:1){
    #block1_traits <- NULL
    block2_traits <- NULL
	
    if(s2 == strip_number){
        # read data for set 1 (first block)
        block1_traits <- strip_list[[s2]]
        cat("Reading (and eventually harmon) ALL TRAITS of set 1 (block", s2, "), has traits:",
        paste(block1_traits, collapse=", "), "\n")

        for (trait in block1_traits) {
          # I should probably make this part a function
	  results[[trait]] <- process_trait(trait, snps_in_ref,"beta_hat","se","af",gwas_info)  
        }
    }

str(results)
stop('you told me to stop here')

    if(s2 == strip_number + 1){
        # read set 2 (second block)
        block2_traits <- strip_list[[s2]]
        cat("Reading (and eventually harmon) ALL TRAITS of set 2 (block", s2, "), has traits:",
        paste(block2_traits, collapse=", "), "\n")
        
        for (trait in block2_traits) {
          results[[trait]] <- process_trait(trait, snps_in_ref)          
        }
    
    
    } else if(s2 > strip_number + 1){
        # drop old set 2 (previous block), read new set 2 (current block)
        prev_block2_traits <- strip_list[[s2 - 1]]
        block2_traits <- strip_list[[s2]]
        cat("Dropping ALL TRAITS of previous set 2 (block", s2-1, "), has traits:",
        paste(prev_block2_traits, collapse=", "), "\n")
        
        cat("Reading (and eventually harmon) ALL TRAITS of new set 2 (block", s2, "), has traits:",
        paste(block2_traits, collapse=", "), "\n")
    
        for (trait in block2_traits) {
          results[[trait]] <- process_trait(trait, snps_in_ref)
        }
    }
    
    str(results)

    # Add within-block comparisons for block1_traits
    #if (!is.null(block1_traits) && is.null(block2_traits)) {
    # need to add handling for blocks which have traits w/ no snps left
    #    within_comparisons <- t(combn(block1_traits, 2))
    #    colnames(within_comparisons) <- c("trait1", "trait2")
    #    within_comparisons <- as.data.frame(within_comparisons, stringsAsFactors = FALSE)
    #    cat("Unique within-block comparisons for block", s2, ":\n")
    #    print(within_comparisons)
        #results[[paste0("block_", s2, "_within")]] <- within_comparisons
	# need N, median of samples sizes, and Zs
    #    l2 <- as.numeric(scan(pipe("awk -F'\t' 'NR>1{print $6}' chr1_8traits_snps_in_ld_file.tsv"), what="character"))
    #    print('l2:')
#	print(head(l2))
        
	#`print(head(Z_hat))
#	result <- R_ldsc(
#           Z_hat = Z_hat,
#           ldscores = l2,
#           N = N,
#           make_well_conditioned = FALSE,
#           comparisons = comparisons
#`:        )
#    }

    ### Pairwise comparison if both blocks are defined:
    if(!is.null(block1_traits) && !is.null(block2_traits)) {
        comparisons <- expand.grid(trait1 = block1_traits,
                                   trait2 = block2_traits,
                                   stringsAsFactors = FALSE)
        cat("Pairwise comparisons between block", strip_number, "and", s2, ":\n")
	print(comparisons)

        # Get unique trait codes
        all_traits <- unique(c(block1_traits, block2_traits))
        # Map each trait name to its index in all_traits
        num_trait1 <- match(comparisons$trait1, all_traits)
        num_trait2 <- match(comparisons$trait2, all_traits)
        # Build numeric data frame (or matrix)
        num_compar <- data.frame(trait1 = num_trait1,
                         trait2 = num_trait2)
	print(num_compar)

        # need N, median of samples sizes, and Zs
        l2 <- as.numeric(scan(pipe(sprintf("awk -F'\t' 'NR>1{print $2}' %s", snps_in_ref)), what="character"))
        print('l2:')
        print(head(l2))
        
        med_ss_vec <- sapply(results, function(x) x$med_ss)
        Z_hat <- as.data.frame(lapply(results, function(x) x$Z))	

        result <- R_ldsc(
           Z_hat = Z_hat,
           ldscores = l2,
           N = med_ss_vec,
           ld_size = M,
	   make_well_conditioned = FALSE,
           comparisons = num_compar
        )
	print(result)

	# Store each result in your ldsc_results list
        block_comp_name <- paste0("block_", strip_number, "_vs_block_", s2)
        ldsc_results[[block_comp_name]] <- list("ldsc"=result,"str_compare"=comparisons,"num_compare"=num_compar)
    }
}

str(results)
str(ldsc_results)

saveRDS(ldsc_results, file=ldsc_out)



