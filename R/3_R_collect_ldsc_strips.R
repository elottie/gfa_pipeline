# snakemake rule
#nstrips=2
#rule R_ldsc_collect:
#    input: ldsc_strip_res = expand(data_dir + "{{prefix}}_R_estimate.R_ldsc.{strip_num}.RDS", strip_num = range(1, nstrips))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    script: "R/collect_ldsc_strips.R"

# let it be a list
ldsc_strip_files <- snakemake@input[["ldsc_strip_res"]]
out <- snakemake@output[["out"]]

#ldsc_strip_res <- readRDS("../gfa_data/First8_Mets_ldsc_results.RDS")
    
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
make_big_symm_matrix <- function(inp_list, row_name, col_name, value_name, flatten=TRUE) {

  if (flatten){
    inp_list <- do.call(c, ldsc_strip_res)
  }

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

big_R_ldsc <- make_big_symm_matrix(ldsc_strip_res, row_name = "trait1", col_name = "trait2", value_name = "intercept")

# here add projection to positive definite

print(big_R_ldsc)

saveRDS(big_R_ldsc, file=out)
