suppressPackageStartupMessages({
  library(data.table)
  library(GFA)
})

source('R/manhattan_helpers.R')

# ============================================================
# Snakemake inputs, outputs, wildcards, and parameters
# ============================================================

gwas_info_file <- snakemake@input[["gwas_info"]]  # make sure this is uncorr_info_input
snp_files <- as.character(snakemake@input[["snp_list"]])
gfa_obj <- snakemake@input[["gfa_obj"]]

trait <- snakemake@wildcards[["trait"]]
analyses <- snakemake@params[["snp_sets"]]
stopifnot(all(unlist(analyses) %in% c('snps_all','snps_eligible_for_gfa','snps_used_in_gfa')))

#output_plot <- snakemake@output[["plot"]]
#output_peaks <- snakemake@output[["peaks"]] 

plot_files <- setNames(
  snakemake@output[["plots"]],
  analyses
)
peak_files <- setNames(
  snakemake@output[["peaks"]],
  analyses
)

genome_build <- as.character(snakemake@params[["genome_build"]])
significance_threshold <- as.numeric(snakemake@params[["sig_thresh_neg_log_10_p"]])
peak_window <- as.numeric(snakemake@params[["peak_window"]])

message("Trait: ", trait)
message("Analyses to plot: ", paste(unlist(analyses), collapse = ", "))
message("Genome build: ", genome_build)
message("Significance threshold: ", significance_threshold)
message("Peak window: ", peak_window, " bp")

workdir <- paste0("7_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Read gwas_info_file and select this trait
# ============================================================

gwas_info <- fread(gwas_info_file)

trait_row <- gwas_info[gwas_info[["name"]] == trait]

# ============================================================
# Get the trait-specific file and column names
# ============================================================

trait_file <- as.character(trait_row[["raw_data_path"]][1])

snp_column <- as.character(trait_row[["snp"]][1])
chrom_column <- as.character(trait_row[["chrom"]][1])
position_column <- as.character(trait_row[["pos"]][1])
# just gather beta and se and compute log10 pvals myself
beta_column <- as.character(trait_row[["beta_hat"]][1])
se_column <- as.character(trait_row[["se"]][1])
# I really don't need A1 and A2, but they are required to use gwas_format
A1_column <- as.character(trait_row[["A1"]][1])
A2_column <- as.character(trait_row[["A2"]][1])

message("Trait file: ", trait_file)
message("SNP column: ", snp_column)
message("Chromosome column: ", chrom_column)
message("Position column: ", position_column)
message("Beta column: ", beta_column)
message("Se column: ", se_column)
message("Nuisance A1 column: ", A1_column)
message("Nuisance A2 column: ", A2_column)

# ============================================================
# Check and read the required trait_file columns
# ============================================================

required_trait_file_columns <- unique(c(
  snp_column,
  chrom_column,
  position_column,
  beta_column,
  se_column,
  A1_column,
  A2_column
))

trait_reader <- if (grepl("\\.gz$", trait_file, ignore.case = TRUE)) {
  sprintf("zcat %s", shQuote(trait_file))
} else {
  sprintf("cat %s", shQuote(trait_file))
}

extract_snps_from_trait <- sprintf(
  'BEGIN { FS = OFS = "\\t" }
   NR == FNR {
       extract[$1] = 1
       next
   }
   FNR == 1 {
       for (i = 1; i <= NF; i++) {
           if ($i == "%s") snp_col = i
       }
       if (!snp_col) {
           print "SNP column not found" > "/dev/stderr"
           exit 1
       }
       print
       next
   }
   $(snp_col) in extract {
       print
   }',
  snp_column
)

if ("snps_all" %in% analyses){

  all_dt <- fread(
    cmd = trait_reader,
    select = required_trait_file_columns
  )

  print('head of all_dt after reading in:')
  print(head(all_dt))

  all_dt <- gwas_format(all_dt,
			snp=snp_column, beta_hat=beta_column, se=se_column, A1=A1_column, A2=A2_column,
			chrom=chrom_column, pos=position_column, compute_pval=FALSE)

  print('head of all_dt after cleaning w/ gwas format:')
  colnames(all_dt) <- c(chrom_column,position_column,snp_column,A1_column,A2_column,beta_column,se_column,'p','ss','af','row_id')
  print(head(all_dt))
}
if ("snps_eligible_for_gfa" %in% analyses){
  # temp output file defs
  elig_snp_file <- file.path(workdir, "elig_snps_all_chr.txt")

  make_elig_snp_file <- sprintf(
    "awk 'FNR > 1 {print $1}' %s > %s",
    paste(shQuote(snp_files), collapse = " "),
    shQuote(elig_snp_file)
  )

  system(make_elig_snp_file)

  trait_elig_read <- sprintf(
    "%s | awk %s %s -",
    trait_reader,
    shQuote(extract_snps_from_trait),
    shQuote(elig_snp_file)
  )

  elig_dt <- fread(
    cmd = trait_elig_read,
    select = required_trait_file_columns
  )

  print('head of elig_dt after reading in:')
  head(elig_dt)
}
if ("snps_used_in_gfa" %in% analyses){
  # temp output file defs
  used_snp_file <- file.path(workdir, "used_snps_all_chr.txt")

  gfa <- readRDS(gfa_obj)
  used_snps <- gfa$snps
  rm(gfa)  # large object.  we don't need anymore
  gc()

  fwrite(
    data.table(snp = used_snps),
      file = used_snp_file,
      col.names = FALSE,
      quote = FALSE
  )

  trait_used_read <- sprintf(
    "%s | awk %s %s -",
    trait_reader,
    shQuote(extract_snps_from_trait),
    shQuote(used_snp_file)
  )

  used_dt <- fread(
    cmd = trait_used_read,
    select = required_trait_file_columns
  )

  print('head of used_dt after reading in:')
  head(used_dt)
}

# make list of the dts we have and do further work on each of them.  list of datatables does not make whole copy but pointer
analysis_dts <- list()

if (exists("all_dt")) {
  analysis_dts$snps_all <- all_dt
}
if (exists("elig_dt")) {
  analysis_dts$snps_eligible_for_gfa <- elig_dt
}
if (exists("used_dt")) {
  analysis_dts$snps_used_in_gfa <- used_dt
}

# ============================================================
# begin work on each dataset as an analysis
# ============================================================

for (analysis in analyses) {

  message('working on analysis: ',analysis)
  
  dt = analysis_dts[[analysis]]

  # get only the columns needed for the plot. harmonization doesn't matter, we just want p
  dt <- dt[, {
    beta <- as.numeric(get(beta_column))
    se   <- as.numeric(get(se_column))

    valid <- !is.na(beta) & !is.na(se) & se > 0
    z <- beta[valid] / se[valid]

    # Numerically stable calculation of -log10(2 * pnorm(-abs(z)))
    negLog10p <- rep(NA_real_, length(beta))
  
    negLog10p[valid] <- -(log(2) + pnorm(
      -abs(z),
      log.p = TRUE
    )) / log(10)

    .(
      snp   = as.character(get(snp_column)),
      chrom = as.numeric(get(chrom_column)),
      pos   = as.numeric(get(position_column)),
      negLog10p  = negLog10p
    )
  }]
  
  print('head of dt after making ready to plot:')
  print(head(dt))

  # make plot!
  plot_manhattan(
    dt = dt,
    analysis_name = analysis,
    genome_build = genome_build,
    significance_threshold = significance_threshold,
    peak_window = peak_window,
    plot_file = plot_files[[analysis]],
    peaks_file = peak_files[[analysis]]
  )

  gc()

}
# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
