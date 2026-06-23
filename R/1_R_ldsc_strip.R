library(dplyr)
#library(tidyr)
library(purrr)
library(readr)
library(GFA)
# library(stringr)

# added
library(data.table)
library(ps)
#library(bench)
library(jsonlite)

# --- for once I graduate to snakemake ---
#rule R_ldsc_strip:
#    input: snp_list = expand(data_dir + "snp_lists/" + "{{prefix}}_snps_chr{chrom}.tsv", chrom = range(1, 23)),
           #raw_data_input, 
#           gwas_info = info_input,
#           strip_list =  data_dir + "{prefix}_ldsc_strip_list.RDS",
#           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
#           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.{strip_num}.RDS"
#    resources: mem_mb = mem_limit_gb # could adjust resources
#    script: "R/3_R_ldsc_strip.R"


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
#snp_files <- unlist(snakemake@input[["snp_list"]])
#gwas_info <- fread(snakemake@input[["gwas_info"]])
#strip_list <- read_json(snakemake@input[["strip_list"]],simplifyVector=TRUE,simplifyMatrix=FALSE)
#strip_num <- as.numeric(snakemake@wildcards[["strip_num"]])
#out <- snakemake@output[["out"]]
#ld_files <- unlist(snakemake@input[["l2"]])
#m_files <- unlist(snakemake@input[["m"]])
# NEEDS NO MORE L2DIR

snp_files <- sprintf("../gfa_data/snp_lists/First8SnakemakeTest_snps_chr%d.tsv", 1:22)
gwas_info <- fread("../First8_Mets_ForLDSCStrip.csv")
strip_list <- read_json("../gfa_data/First8SnakemakeTest_ldsc_strip_list.json",simplifyVector=TRUE,simplifyMatrix=FALSE)  # list of character vectors
strip_num <- 1
# these go away with snakemake input --
l2_dir <- "/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/"
chroms <- 1:22
# --
m_files <- paste0(l2_dir, chroms, ".l2.M_5_50")
ld_files <- paste0(l2_dir, chroms, ".l2.ldscore.gz")
out <- "../gfa_data/First8SnakemakeTest_ldsc_results.RDS"

# --- source helpful funcs ---
# make awks, sorts, and joins consistent across users
Sys.setenv(LC_ALL = "C")

harmon_helper_path <- "harmon_helpers.R"
ld_ref_helper_path <- "ld_ref_helpers.R"
# eventually need to switch to this
#helper_path <- "R/harmon_helpers.R"
source(harmon_helper_path)
source(ld_ref_helper_path)

curr_ram <- function(label = "") {
  rss <- ps::ps_memory_info(ps::ps_handle())[["rss"]] / 1024^3
  cat(sprintf("%-30s %8.3f GB\n", label, rss))
  invisible(rss)
}

curr_ram('after reading input')

# --- temp workdir for testing cleanliness --
workdir <- paste0("1_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
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
curr_ram('after ld fef concatenation')

# --- other random input setup ---
# if M is num of variants used to compute ld scores, should be constant?
M <- purrr:::map(1, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

# read in gwas_info for trait names
traits <- gwas_info$name

# --- set up blocks for input to ldsc_strip ---
print(paste("Received smart-selected nblocks:", length(strip_list)))
print('These blocks are:')
print(strip_list)

curr_ram('after other random input setup')

# --- strip processing! ---
# read in the ref snps for rownames
# I have already deduplicated them
snps_in_ref <- fread(snps_in_ref_file, header = TRUE, select = 1)[[1]]
n_snps <- length(snps_in_ref)
print(n_snps)
n_traits <- length(traits)

# max block2 size for making size of Z
max_block2_size <- if (length(strip_list) >= 2) max(lengths(strip_list)[-1]) else 0
block1_traits <- strip_list[[strip_num]]
Z_work <- matrix(NA_real_, n_snps, length(block1_traits) + max_block2_size,
                 dimnames = list(snps_in_ref, c(block1_traits, rep("", max_block2_size))))
ss_work <- matrix(NA_real_, n_snps, length(block1_traits) + max_block2_size,
                 dimnames = list(snps_in_ref, c(block1_traits, rep("", max_block2_size))))

l2 <- as.numeric(scan(pipe(sprintf("awk 'NR > 1 {print $2}' %s", snps_in_ref_file)), what="character"))
#print('l2:')
#print(head(l2))
#print(length(l2))
#print(nrow(Z_work))

curr_ram('after strip setup')

ldsc_results <- list()
for(s2 in strip_num:(length(strip_list))){
#for(s2 in strip_num:2){
#for(s2 in strip_num:1){

    if(s2 == strip_num){
        # read data for set 1 (first block)
        #block1_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 1 (block", s2, "), has traits:",
        paste(block1_traits, collapse=", "), "\n")

        for (trait in block1_traits) {
          harmon <- harmon_dat(gwas_info, trait, snps_in_ref_file, return_ss=TRUE)
          
          if (identical(harmon$snps,rownames(Z_work))){
            Z_work[, trait] <- harmon$Z
            ss_work[, trait] <- harmon$ss
          } else {
            stop('rowname snps used in harmon_dat are not the same as rownames of destination Z_work matrix')
          }
          gc()	  
        }
        
        # for ldsc:  add within-block comparisons for block1_traits
        comp_idx <- which(upper.tri(matrix(TRUE, length(block1_traits), length(block1_traits)), diag = TRUE), arr.ind = TRUE)
        comparisons <- data.frame(
          trait1 = block1_traits[comp_idx[, 1]],
          trait2 = block1_traits[comp_idx[, 2]],
          stringsAsFactors = FALSE
        )
        cat("Unique within-block comparisons for block", s2, ":\n")
        print(comparisons)

        curr_ram('after within-block comparison')
    }

    if(s2 > strip_num){
        # read set 2 (second block)
        block2_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 2 (block", s2, "), has traits:",
        paste(block2_traits, collapse=", "), "\n")

        for (j in seq_along(block2_traits)) {
            trait <- block2_traits[j]
            harmon <- harmon_dat(gwas_info, trait, snps_in_ref_file, return_ss=TRUE)
            
	    if (identical(harmon$snps,rownames(Z_work))){
              Z_work[, length(block1_traits) + j] <- harmon$Z
              ss_work[, length(block1_traits) + j] <- harmon$ss
              colnames(Z_work)[length(block1_traits) + j] <- trait
              colnames(ss_work)[length(block1_traits) + j] <- trait
            } else {
              stop('rowname snps used in harmon_dat are not the same as rownames of destination Z_work matrix')
            }
	    gc()   
	} 

        # there may be some columns of Z_work and ss_work that are NA because of the way we defined the matrix to be > size of max block2
        # this is okay.  ldsc_rg will only work with traits that are referenced in comparisons.  NA traits will not be referenced
        
        # for ldsc, pairwise comparison if both blocks are defined:
        comparisons <- expand.grid(trait1 = block1_traits,
                                   trait2 = block2_traits,
                                   KEEP.OUT.ATTRS = FALSE,
                                   stringsAsFactors = FALSE)
        cat("Pairwise comparisons between block", strip_num, "and", s2, ":\n")
        print(comparisons)

        curr_ram('after between-block comparison')
        
        # second to last-strip has special handling to get last triangle for free
        if(strip_num == length(strip_list)-1){
            # read data for set 1 (first block)
            #block1_traits <- strip_list[[s2]]
            cat("Since we are at second-to-last strip, doing last strip, i.e. adding within-block comparisons for very last block ", s2, ", has traits:",
            paste(block2_traits, collapse=", "), "\n")

            # for ldsc:  add within-block comparisons for block1_traits
            comp_idx <- which(upper.tri(matrix(TRUE, length(block1_traits), length(block1_traits)), diag = TRUE), arr.ind = TRUE)
            last_comparisons <- data.frame(
              trait1 = block2_traits[comp_idx[, 1]],
              trait2 = block2_traits[comp_idx[, 2]],
              stringsAsFactors = FALSE
            )
            cat("Unique within-block comparisons for block", s2, ":\n")
            print(last_comparisons)

            # add to earlier comparisons
            comparisons <- rbind(comparisons, last_comparisons)
            comparisons <- unique(comparisons)  # optional, removes exact duplicate rows

            curr_ram('after last within-block comparison')
        }
    }
    
    # need Zs and sample sizes
    print('head of Z_work:')
    print(head(Z_work))
    print(paste('sum of NAs in Z_work:',sum(is.na(Z_work))))
    print('head of ss_work:')
    print(head(ss_work))

    ldsc_result <- R_ldsc(
        Z_hat = Z_work,
        ldscores = l2,
        N = ss_work,
        ld_size = M,
        make_well_conditioned = FALSE,
        comparisons = comparisons
    )
    print(ldsc_result)

    curr_ram('after ldsc_result')

    # Store each result in your ldsc_results list
    block_comp_name <- paste0("block_", strip_num, "_vs_block_", s2)
    ldsc_results[[block_comp_name]] <- ldsc_result

    curr_ram('after storing ldsc result')

    gc()
}

str(ldsc_results)
saveRDS(ldsc_results, file=out)

curr_ram('after saving results')

# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
