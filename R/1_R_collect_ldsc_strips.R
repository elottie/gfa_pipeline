library(data.table)

# snakemake rule
#nstrips=2
#rule R_ldsc_collect:
#    input: ldsc_strip_res = expand(data_dir + "{{prefix}}_R_estimate.R_ldsc.{strip_num}.RDS", strip_num = range(1, nstrips))
#    output: out = data_dir + "{prefix}_R_estimate.R_ldsc.RDS"
#    script: "R/collect_ldsc_strips.R"

# let it be a list
gwas_info_file <- snakemake@input[["gwas_info"]]
ldsc_strip_files <- snakemake@input[["ldsc_strip_res"]]
cor_cutoff <- snakemake@params[["cor_cutoff"]]
cond_num <- snakemake@params[["cond_num"]]
out <- snakemake@output[["out"]]
uncorr_info <- snakemake@output[["uncorr_info"]]

#gwas_info_file <- "../First8_Mets_ForLDSCStrip.csv"
#ldsc_strip_files <- "../gfa_data/First8SnakemakeTest_ldsc_results.RDS"
#cor_cutoff <- 0.97
#cond_num <- 100
#out <- "../gfa_data/First8SnakemakeTest_collect_R_strips.RDS"
#uncorr_info <- file.path(dirname(out), sub("\\.csv$", "_uncorr_traits.csv", gwas_info_file))

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

make_big_symm_matrix <- function(inp_list, row_name, col_name, value_name, flatten=TRUE, cor_cutoff=1) {

  if (flatten){
    inp_list <- unlist(inp_list, recursive = FALSE, use.names = TRUE)
  }

  # 1. Combine all pairwise correlation/intercept results into one table
  all_pair_table <- rbindlist(lapply(inp_list, function(item) {
    data.table(
      trait_1 = as.character(item[[row_name]]),
      trait_2 = as.character(item[[col_name]]),
      pair_val = as.numeric(item[[value_name]])
    )
  }), fill = TRUE)

  # 2. Extract high-correlation pairs
  # mess it up
  print(all_pair_table)
#  all_pair_table[c(2,4),'pair_val'] <- 0.98
#  all_pair_table[c(7,8),'pair_val'] <- -0.98
#  print(all_pair_table)
  
  high_pairs <- all_pair_table[
    !is.na(pair_val) & trait_1 != trait_2 & abs(pair_val) > cor_cutoff,
    .(
      # if have C11 and C12, always represent as C11-C12, never C12-C11
      trait_1 = pmin(trait_1, trait_2),
      trait_2 = pmax(trait_1, trait_2),
      pair_val,
      abs_pair_val = abs(pair_val)
    )
  ][
    # sort rows by descending absolute value.  if dup pairs, keep the row with highest abs_pair_val (though here should be identical)
    order(-abs_pair_val),
    .SD[1],
    by = .(trait_1, trait_2)
  ] 
  
  print(high_pairs)

  # drop traits w/ most corrs first, descend
  drop_traits <- character(0)

  drop_history <- data.table(
    step = integer(),
    trait = character(),
    n_high_corr = integer(),
    max_pair_val = numeric(),
    mean_pair_val = numeric(),
    correlated_traits = character()
  )

  remaining_high_pairs <- copy(high_pairs)

  step <- 1L

  while (nrow(remaining_high_pairs) > 0) {

    # Make long table: one row per trait-pair direction
    cor_edges_long <- rbind(
      remaining_high_pairs[, .(trait = trait_1, correlated_trait = trait_2, pair_val = pair_val)],
      remaining_high_pairs[, .(trait = trait_2, correlated_trait = trait_1, pair_val = pair_val)]
    )

    # Count high-correlation links among currently remaining traits
    cor_counts <- cor_edges_long[
      ,
      .(
	n_high_corr = .N,
	max_pair_val = pair_val[which.max(abs(pair_val))],
	mean_pair_val = mean(abs(pair_val), na.rm = TRUE),
	correlated_traits = paste(sort(unique(correlated_trait)), collapse = ", ")
      ),
      by = trait
    ]


    # Pick trait to drop:
    # 1. highest number of high-correlation links
    # 2. if tied, highest absolute mean correlation
    # 3. if tied, alphabetical trait name for reproducibility
    trait_to_drop_row <- cor_counts[
      order(-n_high_corr, -abs(mean_pair_val), trait)
    ][1]

    trait_to_drop <- trait_to_drop_row$trait

    # Save history, now including correlated traits
    drop_history <- rbind(
      drop_history,
      data.table(
        step = step,
        trait = trait_to_drop,
        n_high_corr = trait_to_drop_row$n_high_corr,
        max_pair_val = trait_to_drop_row$max_pair_val,
        mean_pair_val = trait_to_drop_row$mean_pair_val,
        correlated_traits = trait_to_drop_row$correlated_traits
      )
    )

    drop_traits <- c(drop_traits, trait_to_drop)

    # Remove all pairs involving dropped trait
    remaining_high_pairs <- remaining_high_pairs[
      trait_1 != trait_to_drop & trait_2 != trait_to_drop
    ]

    step <- step + 1L
  }

  drop_traits <- unique(drop_traits)

  print(drop_history)
  print(drop_traits)

  # Remove all rows involving dropped traits
  uncorr_pair_table <- all_pair_table[!(trait_1 %chin% drop_traits | trait_2 %chin% drop_traits)]

  # 1. Get all kept unique trait names
  uncorr_traits <- unique(c(uncorr_pair_table$trait_1, uncorr_pair_table$trait_2))
  n_traits <- length(uncorr_traits)

  # 3. Preallocate final matrix once
  out <- matrix(
    NA_real_,
    nrow = n_traits,
    ncol = n_traits,
    dimnames = list(uncorr_traits, uncorr_traits)
  )

  # 4. fill matrix
  i <- match(uncorr_pair_table$trait_1, uncorr_traits)
  j <- match(uncorr_pair_table$trait_2, uncorr_traits)

  # cbind so it does each i,j pair, not looks for a block
  out[cbind(i, j)] <- uncorr_pair_table$pair_val
  out[cbind(j, i)] <- uncorr_pair_table$pair_val

 return(list(symm_mat=out, drop_history=drop_history))
}

symm_mat_res <- make_big_symm_matrix(ldsc_strip_res, row_name = "trait1", col_name = "trait2", value_name = "intercept", cor_cutoff=cor_cutoff)
big_ldsc_se <- symm_mat_res$symm_mat
drop_traits <- symm_mat_res$drop_history
print(big_ldsc_se)

# here add projection to positive definite
pos_def_se <- Matrix::nearPD(big_ldsc_se, corr = FALSE, keepDiag = TRUE,
                         posd.tol = 1/cond_num)$mat
pos_def_se <- as.matrix(pos_def_se)
print(pos_def_se)

# now we have subset of traits we want to use in later analysis.  the easiest thing to do since later analysis reads gwas_info is to write new gwas_info
stopifnot(identical(rownames(big_ldsc_se),colnames(big_ldsc_se)))
uncorr_traits <- colnames(big_ldsc_se)
gwas_info_uncorr <- gwas_info[name %in% uncorr_traits]

# --- save outputs ---
saveRDS(pos_def_se, file=out)
fwrite(gwas_info_uncorr, uncorr_info)
fwrite(drop_traits, sub(
  "_uncorr_traits\\.csv$",
  "_dropped_corr_traits.tsv",
  uncorr_info, sep='\t'
))
