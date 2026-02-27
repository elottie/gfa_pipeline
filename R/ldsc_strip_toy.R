# library(dplyr)
# library(purrr)
# library(readr)
library(GFA)
# library(stringr)

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
Z_hat <- matrix(c(
  1.2, 0.8,   # SNP1
  1.1, 0.7,   # SNP2
  0.9, 0.5,   # SNP3
  1.0, 0.6    # SNP4
), nrow = 4, byrow = TRUE)
# LD scores (numeric vector) for each SNP.  come from LD ref panel, not me
ldscores <- c(1.5, 2.0, 1.8, 2.2)
# LD size, num of variants used to compute LD scores, right now arbitrary scalar
ld_size <- 4
# N: sample sizes per trait (if equal across SNPs, can be a vector)
N <- c(1000, 1200)
blocks <- NULL
ncores <- 1
# the new option, comparisons of traits
comparisons <- matrix(c(1,2), ncol=2)
print(comparisons)
# tell Jean abt 2-col err!

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

# --- set up blocks for input to ldsc_strip.  one day, a rule makes? ---
get_set <- function(mat, block_size, i, j, trait_names) {
  row_start <- (i - 1) * block_size + 1
  row_end   <- min(i * block_size, nrow(mat))
  col_start <- (j - 1) * block_size + 1
  col_end   <- min(j * block_size, ncol(mat))
  block <- mat[row_start:row_end, col_start:col_end]
  row_names <- trait_names[row_start:row_end]
  col_names <- trait_names[col_start:col_end]
  list(block=block, row_names=row_names, col_names=col_names)
}
#get_triangle <- function(mat, block_size, i, trait_names) {
#  set <- get_set(mat, block_size, i, i, trait_names)
#  mask <- upper.tri(set$block)
#  inds <- which(mask, arr.ind=TRUE)
#  data.frame(
#    trait1 = set$row_names[inds[,1]],
#    trait2 = set$col_names[inds[,2]],
#    value  = set$block[mask]
#  )
#}

# example traits, expand them out to a matrix for me to understand what is happening
trait_names <- paste0("Trait", 1:8)
mat <- matrix(NA, ncol=8, nrow = 8)
dimnames(mat) <- list(trait_names, trait_names)

block_size <- 2
n_blocks <- ceiling(length(trait_names) / block_size)   # 4 in this example

all_blocks <- list()
# make the blocks and store them in all_blocks.  treat triangles like blocks for now
for (i in 1:n_blocks) {
  for (j in 1:n_blocks) {
    if (j > i) { # Above diagonal blocks: full pairs
      block_name <- paste0("block_", i, "_", j)
      all_blocks[[block_name]] <- get_set(mat, block_size, i, j, trait_names)
    } else if (j == i) { # On the diagonal: only store upper triangle pairs
      block_name <- paste0("triangle_", i)
      #all_blocks[[block_name]] <- get_triangle(mat, block_size, i, trait_names)
      all_blocks[[block_name]] <- get_set(mat, block_size, i, j, trait_names)
    }
    # (j < i) is not stored, as those are redundant
  }
}

# show block names and the blocks themselves
names(all_blocks)
all_blocks

# --- create the relevant strip out of all blocks ---

# do I expect the user to provide the right blocks, or do I need to select strip myself?
strip_list <- all_blocks[1:4]
strip_number <- 1

nblocks <- length(strip_list)

# --- strip processing! ---
results <- list()
for(s2 in strip_number:nblocks){
    if(s2 == strip_number){
        # read data for set 1 (first block)
        block_name <- names(strip_list)[s2]
        block1_traits <- strip_list[[s2]]
        is_triangle <- grepl("^triangle_", block_name)

        if(is_triangle) {
            cat("Reading TRIANGLE of set 1 (", block_name, "), has row traits:",
                paste(block1_traits$row_names, collapse=", "),
                "; column traits:",
                paste(block1_traits$col_names, collapse=", "), "\n")
        } else {
            cat("Reading ALL TRAITS of set 1 (", block_name, "), has row traits:",
                paste(block1_traits$row_names, collapse=", "),
                "; column traits:",
                paste(block1_traits$col_names, collapse=", "), "\n")
        }
    
    }
    if(s2 == strip_number + 1){
        # read set 2 (second block)
        block2_traits <- strip_list[[s2]]
        cat("Reading JUST COLUMN TRAITS of set 1 (block", s2, "), has row traits:",
        paste(block2_traits$row_names, collapse=", "),
        "; column traits:",
        paste(block2_traits$col_names, collapse=", "), "\n")
        
    } else if(s2 > strip_number + 1){
        # drop old set 2 (previous block), read new set 2 (current block)
        prev_block2_traits <- strip_list[[s2 - 1]]
        new_block2_traits <- strip_list[[s2]]
        cat("Dropping JUST COLUMN TRAITS of previous set 2 (block", s2-1, "), has row traits:",
        paste(prev_block2_traits$row_names, collapse=", "),
        "; column traits:",
        paste(prev_block2_traits$col_names, collapse=", "), "\n")
        
        cat("Reading JUST COLUMN TRAITS of new set 2 (block", s2, "), has row traits:",
        paste(new_block2_traits$row_names, collapse=", "),
        "; column traits:",
        paste(new_block2_traits$col_names, collapse=", "), "\n")
    }
    
    # Example calculation
    #results[[s2]] <- list(block_num = s2, snps = strip_list[[s2]], n_snps = length(strip_list[[s2]]))
}

#print(results)

#saveRDS(ret, file=out)



