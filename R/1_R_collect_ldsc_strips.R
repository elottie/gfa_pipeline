library(data.table)

# snakemake rule
#nstrips=2
#rule R_ldsc_collect:
#    input: ldsc_strip_res = expand(data_dir + "{{prefix}}_R_estimate.R_ldsc.{strip_num}.RDS", strip_num = range(1, nstrips))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    script: "R/collect_ldsc_strips.R"

# let it be a list
#gwas_info_file <- snakemake@input[["gwas_info"]]
#ldsc_strip_files <- snakemake@input[["ldsc_strip_res"]]
#cor_cutoff <- snakemake@params[["cor_cutoff"]]
#cond_num <- snakemake@params[["cond_num"]]
#out <- snakemake@output[["out"]]
#uncorr_info <- snakemake@output[["uncorr_info"]]

#gwas_info <- fread(gwas_info_file)

gwas_info_file <- "../First8_Mets_ForLDSCStrip.csv"
ldsc_strip_files <- "../gfa_data/First8SnakemakeTest_ldsc_results.RDS"
cor_cutoff <- 0.97
cond_num <- 100
out <- "../gfa_data/First8SnakemakeTest_collect_R_strips.RDS"
uncorr_info <- sub("\\.csv$", "_uncorr_traits.csv", gwas_info_file)

gwas_info <- fread(gwas_info_file)

#GFA:::make_symm_matrix(res, row_name = "trait1", col_name = "trait2", value_name = "intercept")

#List of 2
# $ block_3_vs_block_3:'data.frame':     1 obs. of  3 variables:
#  ..$ trait1   : chr "C999912844"
#  ..$ trait2   : chr "C999913431"
#  ..$ intercept: num 0.062
# $ block_3_vs_block_4:'data.frame':     5 obs. of  3 variables:
#  ..$ trait1   : chr [1:5] "C999912844" "C999913431" "C999912844" "C999913431" ...
#  ..$ trait2   : chr [1:5] "C999913553" "C999913553" "C999916935" "C999916935" ...
#  ..$ intercept: num [1:5] 0.05315 0.12661 -0.00293 -0.01567 0.08051

#for (result_num in 1:length(ldsc_strip_res)){
#  print(result_num)
#  result <- ldsc_strip_res[[result_num]]
#  print(head(result))
#
#  symm_res <- GFA:::make_symm_matrix(result, row_name = "trait1", col_name = "trait2", value_name = "intercept")
#  print(head(symm_res))
#}

ldsc_strip_res <- lapply(ldsc_strip_files, readRDS)

# this is just slightly altered version of make_symm_matrix so we don't have to do that and join pieces
make_big_symm_matrix <- function(inp_list, row_name, col_name, value_name, flatten=TRUE, cor_cutoff=1) {

  if (flatten){
    inp_list <- unlist(inp_list, recursive = FALSE, use.names = TRUE)
  }

  # mess it up
  #print(inp_list)
  inp_list$block_1_vs_block_1[c(1,2),'intercept'] <- 0.98
  #inp_list$block_1_vs_block_2[c(3,6,9),'intercept'] <- 0.98
  #print(inp_list)

  # 0. greedy cor clustering.  if traits a and b have intercept > cor_cutoff, remove b
  drop_traits <- character(0)

  for (item in inp_list) {
    row_item <- as.character(item[[row_name]])
    col_item <- as.character(item[[col_name]])
    pair_val <- as.numeric(item[[value_name]])

    #print(pair_val)

    high_corr <- !is.na(pair_val) & (row_item != col_item) & (pair_val > cor_cutoff)

    if (any(high_corr)) {
      for (i in which(high_corr)) {
        high_corr_item_a <- row_item[i]
        high_corr_item_b <- col_item[i]

        # If neither trait has already been dropped, keep the first trait and drop the second.
        if (!(high_corr_item_a %in% drop_traits) && !(high_corr_item_b %in% drop_traits)) {
          drop_traits <- c(drop_traits, high_corr_item_b)
        }
      }
    }
  }

  drop_traits <- unique(drop_traits)

  #print(drop_traits)

  # Remove all rows involving dropped traits
  if (length(drop_traits) > 0) {
    inp_list <- lapply(inp_list, function(item) {
      keep <- !(
        as.character(item[[row_name]]) %in% drop_traits |
          as.character(item[[col_name]]) %in% drop_traits
      )

      item[keep, , drop = FALSE]
    })

    # Remove now-empty data.frames
    inp_list <- inp_list[vapply(inp_list, nrow, integer(1)) > 0]
  }

  print(inp_list)

  # 1. Get all unique trait names
  cols <- unique(unlist(
    lapply(inp_list, function(x) {
      c(as.character(x[[row_name]]), as.character(x[[col_name]]))
    }),
    use.names = FALSE
  ))

  n_cols <- length(cols)

  # 2. Create lookup: trait name -> matrix index
  col_index <- seq_len(n_cols)
  names(col_index) <- cols

  # 3. Preallocate final matrix once
  out <- matrix(
    NA_real_,
    nrow = n_cols,
    ncol = n_cols,
    dimnames = list(cols, cols)
  )

  # 4. Fill matrix
  for (item in inp_list) {
    ri <- col_index[as.character(item[[row_name]])]
    ci <- col_index[as.character(item[[col_name]])]
    v <- item[[value_name]]

    out[cbind(ri, ci)] <- v
    out[cbind(ci, ri)] <- v
  }

 return(out)
}

big_ldsc_se <- make_big_symm_matrix(ldsc_strip_res, row_name = "trait1", col_name = "trait2", value_name = "intercept", cor_cutoff=cor_cutoff)
print(big_ldsc_se)

# here add projection to positive definite
pos_def_se <- Matrix::nearPD(big_ldsc_se, corr = FALSE, keepDiag = TRUE,
                         posd.tol = 1/cond_num)$mat
print(pos_def_se)

# now we have subset of traits we want to use in later analysis.  the easiest thing to do since later analysis reads gwas_info is to write new gwas_info
stopifnot(identical(rownames(big_ldsc_se),colnames(big_ldsc_se)))
uncorr_traits <- colnames(big_ldsc_se)
gwas_info_uncorr <- gwas_info[name %in% uncorr_traits]

# --- save outputs ---
saveRDS(pos_def_se, file=out)
fwrite(gwas_info_uncorr, uncorr_info)
