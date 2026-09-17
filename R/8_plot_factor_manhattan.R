suppressPackageStartupMessages({
  library(data.table)
  library(GFA)
})

source('R/manhattan_helpers.R')

# ============================================================
# Snakemake inputs, outputs, wildcards, and parameters
# ============================================================

loading_files <- as.character(snakemake@input[["loading_list"]])
gfa_obj <- snakemake@input[["gfa_obj"]]
factor <- as.integer(snakemake@wildcards[["factor"]])

# don't treat factors as wildcard.  could to save time
# this analysis can only include snps that are eligible for gfa or used in gfa
analyses <- snakemake@params[["snp_sets"]]
stopifnot(all(unlist(analyses) %in% c('snps_all','snps_eligible_for_gfa','snps_used_in_gfa')))
analyses_and_factor <- paste0(analyses,"_factor_",factor)

plot_files <- setNames(
  snakemake@output[["plots"]],
  analyses_and_factor
)
peak_files <- setNames(
  snakemake@output[["peaks"]],
  analyses_and_factor
)

genome_build <- as.character(snakemake@params[["genome_build"]])
sig_thresh_negLog10p <- as.numeric(snakemake@params[["sig_thresh_negLog10p"]])
peak_window <- as.numeric(snakemake@params[["peak_window"]])

message("Plotting manhattans for GFA factor: ", factor)
message("Analyses to plot: ", paste(unlist(analyses), collapse = ", "))
message("Genome build: ", genome_build)
message("Significance threshold: ", sig_thresh_negLog10p)
message("Peak window: ", peak_window, " bp")

workdir <- paste0("8_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Check and read the required loading_file columns
# ============================================================

z_column <- paste0('factor',factor,'.z') 
required_loading_file_columns <- unique(c(
  "snp",
  "chrom",
  "pos",
  z_column
))

# there is no handling for snps_all for the factors
elig_dt <- rbindlist(
  lapply(
    loading_files,
    fread,
    select = required_loading_file_columns
  ),
  use.names = TRUE,
  fill = FALSE
)
      
print('head of elig_dt after reading in:')
print(head(elig_dt))

if ("snps_used_in_gfa" %in% analyses){
  gfa <- readRDS(gfa_obj)
  used_snps <- gfa$snps
  rm(gfa)  # large object.  we don't need anymore
  gc()

  used_dt <- elig_dt[snp %chin% used_snps]  
  
  if (!("snps_eligible_for_gfa" %in% analyses)){
    rm(elig_dt)
    gc()
  }

  print('head of used_dt after reading in:')
  print(head(used_dt))
}

# make list of the dts we have and do further work on each of them.  list of datatables does not make whole copy but pointer
analysis_dts <- list()

if (exists("elig_dt")) {
  analysis_dts$snps_eligible_for_gfa <- elig_dt
}
if (exists("used_dt")) {
  analysis_dts$snps_used_in_gfa <- used_dt
}

# ============================================================
# begin work on each dataset as an analysis
# ============================================================

for (analysis in analyses_and_factor) {
  message('working on analysis: ',analysis)

  analysis_no_factor <- gsub(paste0("_factor_",factor),"",analysis)
  dt = analysis_dts[[analysis_no_factor]]

  # get only the columns needed for the plot. harmonization doesn't matter, we just want p
  dt <- dt[, {
    z <- as.numeric(get(z_column))
    valid <- !is.na(z)

    # Numerically stable calculation of -log10(2 * pnorm(-abs(z)))
    negLog10p <- rep(NA_real_, length(z))
  
    negLog10p[valid] <- -(log(2) + pnorm(
      -abs(z[valid]),
      log.p = TRUE
    )) / log(10)

    .(
      snp   = as.character(snp),
      chrom = as.numeric(chrom),
      pos   = as.numeric(pos),
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
    sig_thresh_negLog10p = sig_thresh_negLog10p,
    peak_window = peak_window,
    plot_file = plot_files[[analysis]],
    peaks_file = peak_files[[analysis]]
  )
  gc()
}

# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)
