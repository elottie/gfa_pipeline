library(dplyr)
#library(tidyr)
library(purrr)
library(readr)
library(GFA)
# library(stringr)

# added
library(data.table)
#library(bench)
library(jsonlite)

# --- for once I graduate to snakemake ---
# sample size affects genetic covariance and h2 but not intercept or genetic correlation
#snp_files <- unlist(snakemake@input[["snp_list"]])
#gwas_info <- read_csv(snakemake@input[["gwas_info"]])
#strip_list <- readRDS(snakemake@input[["strip_file"]])
#strip_num <- as.numeric(snakemake@wildcard[["strip_num"]])
#out <- snakemake@output[["out"]]
#ld_files <- unlist(snakemake@input[["l2"]])
#m_files <- unlist(snakemake@input[["m"]])

snp_files <- sprintf("../gfa_data/5e5Sig_Herit_Mets8_snps_chr%d.tsv", 1:22)
gwas_info <- fread("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv")
strip_list <- fromJSON("ldsc_trait_sets.json", simplifyVector = FALSE)  # list of character vectors
strip_num <- 3
# these go away with snakemake input --
l2_dir <- "/nfs/turbo/sph-jvmorr/ld_reference_files/ldsc_reference_files/eur_w_ld_chr/"
chroms <- 1:22
# --
m_files <- paste0(l2_dir, chroms, ".l2.M_5_50")
ld_files <- paste0(l2_dir, chroms, ".l2.ldscore.gz")
out <- "../gfa_data/5e5Sig_Herit_Mets8_ldsc_results.RDS"

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

#rule R_ldsc_strip:
#    input: snp_list = data_dir + "{prefix}_snps_chr{chrom}.tsv",
#           #raw_data_input,
#           strip_list =  data_dir + "{prefix}_trait_sets.json"
#           m = expand(l2_dir + "{chrom}.l2.M_5_50", chrom = range(1, 23)),
#           l2 = expand(l2_dir + "{chrom}.l2.ldscore.gz", chrom = range(1, 23))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.{strip_num}.RDS"
#    script: "R/ldsc_strip_toy.R"

# --- source helpful funcs ---
helper_path <- "harmon_helpers.R"
# eventually need to switch to this
#helper_path <- "R/harmon_helpers.R"
source(helper_path)

# --- temp workdir for testing cleanliness --
workdir <- paste0("3_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# temp output file defs
snps_in_ref <- file.path(workdir, "snps_in_ld_file.tsv")

# --- get just SNPs present in all traits AND ld ref files ---
print('filtering universal_snps.txt to only those found in ld ref files')

# All lines from ld_file.tsv where the second column matches a value in the first column of snps_pass_all_filts.txt.
# write header to output
system(sprintf("echo -e 'snp\tl2' > %s", shQuote(snps_in_ref)))

for (chr in 1:22) {
  snp_file <- snp_files[chr]
  ld_file  <- file.path(l2_dir, paste0(chr, ".l2.ldscore.gz"))

  awk_snps_in_ref <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==FNR{snps[$1]=1; next} snps[$2]{print $2, $6}' %s - >> %s",
    shQuote(ld_file),
    shQuote(snp_file),
    shQuote(snps_in_ref)
  )
  system(awk_snps_in_ref)
}
print('completed ld filtering')


# other random input setup ---
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

# --- strip processing! ---
n_snps <- as.integer(system(paste0("wc -l < ",snps_in_ref), intern = TRUE)) - 1L   # number of SNPs after filtering, -1 for header
print(n_snps)
n_traits <- length(traits)

# max block2 size for making size of Z
max_block2_size <- if (length(strip_list) >= 2) max(lengths(strip_list)[-1]) else NA_integer_
block1_traits <- strip_list[[strip_num]]
Z_work <- matrix(NA_real_, n_snps, length(block1_traits) + max_block2_size,
                 dimnames = list(NULL, c(block1_traits, rep("", max_block2_size))))
ss_work <- matrix(NA_real_, n_snps, length(block1_traits) + max_block2_size,
                 dimnames = list(NULL, c(block1_traits, rep("", max_block2_size))))

l2 <- as.numeric(scan(pipe(sprintf("awk 'NR>1{print $2}' %s", snps_in_ref)), what="character"))
#print('l2:')
#print(head(l2))
#print(length(l2))
#print(nrow(Z_work))

ldsc_results <- list()
for(s2 in strip_num:length(strip_list)){
#for(s2 in strip_num:2){
#for(s2 in strip_num:1){
    #Z_hat_b2 <- NULL
    #ss_b2 <- NULL	

    if(s2 == strip_num){
        # read data for set 1 (first block)
        #block1_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 1 (block", s2, "), has traits:",
        paste(block1_traits, collapse=", "), "\n")

        for (trait in block1_traits) {
	  harmon <- harmon_dat(gwas_info, trait, snps_in_ref, return_ss=TRUE)
		
	  Z_work[, trait] <- harmon$Z
          ss_work[, trait] <- harmon$ss
	}
	
	# for ldsc:  add within-block comparisons for block1_traits
	comparisons <- as.data.frame(t(combn(block1_traits, 2)), stringsAsFactors = FALSE)
        names(comparisons) <- c("trait1", "trait2")
        cat("Unique within-block comparisons for block", s2, ":\n")
        print(comparisons)
        # Get unique trait codes and turn string trait names to numbers
        all_traits <- unique(as.character(block1_traits))
    }

    if(s2 > strip_num){
        # read set 2 (second block)
        block2_traits <- strip_list[[s2]]
        cat("Reading and harmon ALL TRAITS of set 2 (block", s2, "), has traits:",
        paste(block2_traits, collapse=", "), "\n")

        for (j in seq_along(block2_traits)) {
            trait <- block2_traits[j]
            harmon <- harmon_dat(gwas_info, trait, snps_in_ref, return_ss=TRUE)
            Z_work[, length(block1_traits) + j] <- harmon$Z
            ss_work[, length(block1_traits) + j] <- harmon$ss
            colnames(Z_work)[length(block1_traits) + j] <- trait
	    colnames(ss_work)[length(block1_traits) + j] <- trait
        }

	# there may be some columns of Z_work and ss_work that are NA because of the way we defined the matrix to be > size of max block2
	# this is okay.  ldsc_rg will only work with traits that are referenced in comparisons.  NA traits will not be referenced
	
	# for ldsc, pairwise comparison if both blocks are defined:
        comparisons <- expand.grid(trait1 = block1_traits,
                                   trait2 = block2_traits,
                                   stringsAsFactors = FALSE)
        cat("Pairwise comparisons between block", strip_num, "and", s2, ":\n")
        print(comparisons)
        # Get unique trait codes and turn string trait names to numbers
        all_traits <- unique(c(as.character(block1_traits), as.character(block2_traits)))
        
	# second to last-strip has special handling to get last triangle for free
	if(strip_num == length(strip_list)-1){
	    # read data for set 1 (first block)
            #block1_traits <- strip_list[[s2]]
            cat("Adding within-block comparisons for very last block ", s2, "), has traits:",
            paste(block2_traits, collapse=", "), "\n")

            # for ldsc:  add within-block comparisons for block1_traits
            last_comparisons <- as.data.frame(t(combn(block2_traits, 2)), stringsAsFactors = FALSE)
            names(last_comparisons) <- c("trait1", "trait2")
            cat("Unique within-block comparisons for block", s2, ":\n")
            print(last_comparisons)

	    # add to earlier comparisons
	    comparisons <- rbind(comparisons, last_comparisons)
            comparisons <- unique(comparisons)  # optional, removes exact duplicate rows

            # Get unique trait codes and turn string trait names to numbers
            all_traits <- unique(c(as.character(all_traits),as.character(block2_traits)))
	}
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

    # need Zs and sample sizes
    print('head of Z_work:')
    print(head(Z_work))
    print('head of ss_work:')
    print(head(ss_work))

    ldsc_result <- R_ldsc(
        Z_hat = Z_work,
        ldscores = l2,
        N = ss_work,
        ld_size = M,
	make_well_conditioned = FALSE,
        comparisons = num_compar
    )
    print(ldsc_result)

    # Store each result in your ldsc_results list
    block_comp_name <- paste0("block_", strip_num, "_vs_block_", s2)
    ldsc_results[[block_comp_name]] <- list("ldsc"=ldsc_result,"trait_name_map"=trait_map)
}

str(ldsc_results)
saveRDS(ldsc_results, file=out)

# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
