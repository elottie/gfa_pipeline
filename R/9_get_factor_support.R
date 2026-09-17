library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

# --- snakemake inputs ---
factor_files <- snakemake@input[["factor_peaks"]]
trait_files <- snakemake@input[["trait_peaks"]]
trait_names <- snakemake@params[["traits"]]
out_file <- snakemake@output[[1]]

# --- read in and checks ---
factor_snps <- read_tsv(factor_files)

stopifnot(length(trait_files) == length(trait_names))
trait_snps <- map2_dfr(
    trait_files,
    trait_names,
    function(file, trait_name) {
        read_tsv(file, show_col_types = FALSE) %>%
            transmute(
                trait = trait_name,
                snp
            )
    }
)

top_factor_snps <- factor_snps %>%
    filter(!is.na(negLog10p)) %>%
    group_by(peak_id) %>%
    slice_max(
        order_by = negLog10p,
        n = 1,
        with_ties = FALSE
    ) %>%
    transmute(
        peak_id,
        top_snp = snp,
        nearest_gene = nearest_gene
    ) %>%
    ungroup()

peak_summary <- factor_snps %>%
    select(peak_id, snp) %>%
    inner_join(trait_snps, by = "snp") %>%
    group_by(peak_id) %>%
    summarise(
        n_supporting_traits = n_distinct(trait),
        n_supporting_snps = n_distinct(snp),
        supporting_traits = paste(sort(unique(trait)), collapse = ","),
        .groups = "drop"
    )

support_table <- factor_snps %>%
    distinct(peak_id) %>%
    left_join(top_factor_snps, by = "peak_id") %>%
    left_join(peak_summary, by = "peak_id") %>%
    mutate(
        n_supporting_traits = replace_na(n_supporting_traits, 0L),
        n_supporting_snps = replace_na(n_supporting_snps, 0L),
	supporting_traits = replace_na(supporting_traits, "")
    )

write_tsv(support_table,out_file)
