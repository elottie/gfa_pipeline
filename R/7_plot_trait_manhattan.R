suppressPackageStartupMessages({
  library(data.table)
  library(CMplot)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# ============================================================
# Snakemake inputs, outputs, wildcards, and parameters
# ============================================================

gwas_info_file <- snakemake@input[["gwas_info"]]  # make sure this is uncorr_info_input
snp_files <- as.character(snakemake@input[["snp_list"]])
trait <- snakemake@wildcards[["trait"]]

output_plot <- snakemake@output[["plot"]]
output_peaks <- snakemake@output[["peaks"]]  # eventually can get rid of

genome_build <- as.character(
  snakemake@params[["genome_build"]]
)

significance_threshold <- as.numeric(
  snakemake@params[["significance"]]
)

peak_window <- as.numeric(
  snakemake@params[["peak_window"]]
)

message("Trait: ", trait)
message("GWAS info file: ", gwas_info_file)
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

message("Trait file: ", trait_file)
message("SNP column: ", snp_column)
message("Chromosome column: ", chrom_column)
message("Position column: ", position_column)
message("Beta column: ", beta_column)
message("Se column: ", se_column)


# ============================================================
# Check and read the required trait_file columns
# ============================================================

required_trait_file_columns <- unique(c(
  snp_column,
  chrom_column,
  position_column,
  beta_column,
  se_column
))

# make list of snps to extract from trait file
# temp output file defs
valid_snp_file <- file.path(workdir, "valid_snps_all_chr.txt")

cat_command <- sprintf(
  "awk 'FNR > 1 {print $1}' %s > %s",
  paste(shQuote(snp_files), collapse = " "),
  shQuote(valid_snp_file)
)

system(cat_command)

trait_reader <- if (grepl("\\.gz$", trait_file, ignore.case = TRUE)) {
  sprintf("zcat %s", shQuote(trait_file))
} else {
  sprintf("cat %s", shQuote(trait_file))
}

awk_program <- sprintf(
  'BEGIN { FS = OFS = "\\t" }
   NR == FNR {
       valid[$1] = 1
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
   $(snp_col) in valid {
       print
   }',
  snp_column
)

cmd <- sprintf(
  "%s | awk %s %s -",
  trait_reader,
  shQuote(awk_program),
  shQuote(valid_snp_file)
)

trait_dt <- fread(
  cmd = cmd,
  sep = "\t",
  select = required_trait_file_columns

)
print('head of trait_dt after reading in')
head(trait_dt)

# Read only the columns needed for the plot
# here could do selection of only valid snps
#trait_dt <- fread(
#  trait_file,
#  select = required_trait_file_columns
#)

trait_dt <- trait_dt[, {
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
print('head of trait_dt after making nice')
head(trait_dt)

gc()


# ============================================================
# Clean chromosome values
# ============================================================

trait_dt <- trait_dt[chrom %between% c(1, 22)]

# ============================================================
# Clean variant identifiers and filter invalid observations
# ============================================================

# valid snp selection handled by awk earlier

# CMplot uses variant identifiers for highlighting. Make them unique
# in case an rsID occurs more than once.
# rm all snps that ever occur more than once
trait_dt <- trait_dt[, if (.N == 1L) .SD, by = snp]

setorder(trait_dt, chrom, pos)

message("Valid variants: ", nrow(trait_dt))
message(
  "Chromosomes represented: ",
  paste(unique(trait_dt$chrom), collapse = ", ")
)

print('head of trait_dt after removing non uniq')
head(trait_dt)

# ============================================================
# Select significant peaks for gene labeling
# ============================================================
#
# Every valid variant is included in the Manhattan plot.
# peak_window only determines which significant variants receive
# labels. It does not restrict the genomic range being plotted.

candidates <- copy(
  trait_dt[negLog10p >= significance_threshold]
)

# Start with the strongest association
setorder(candidates, -negLog10p)

peaks <- candidates[0]

while (nrow(candidates) > 0L) {

  lead_variant <- candidates[1]

  peaks <- rbind(
    peaks,
    lead_variant,
    use.names = TRUE
  )

  # Remove other significant candidates on the same chromosome
  # that fall within +/- peak_window of this lead variant.
  candidates <- candidates[
    chrom != lead_variant$chrom |
      abs(pos - lead_variant$pos) > peak_window
  ]
}

message("Significant peaks selected for labeling: ", nrow(peaks))


# ============================================================
# Find the nearest gene for each selected peak
# ============================================================

peaks[, `:=`(
  nearest_gene = NA_character_,
  distance_to_gene = NA_integer_
)]

if (nrow(peaks) > 0L) {

  if (genome_build %in% c("GRCh38", "hg38")) {

    suppressPackageStartupMessages({
      library(EnsDb.Hsapiens.v86)
    })

    annotation_database <- EnsDb.Hsapiens.v86

  } else if (genome_build %in% c("GRCh37", "hg19")) {

    suppressPackageStartupMessages({
      library(EnsDb.Hsapiens.v75)
    })

    annotation_database <- EnsDb.Hsapiens.v75

  } else {

    stop(
      "Unsupported genome build: ",
      genome_build,
      ". Expected GRCh38/hg38 or GRCh37/hg19."
    )
  }

  peak_ranges <- GRanges(
    seqnames = peaks$chrom,
    ranges = IRanges(
      start = as.integer(peaks$pos),
      end = as.integer(peaks$pos)
    )
  )

  gene_ranges <- genes(
    annotation_database,
    columns = c("gene_id", "gene_name")
  )

  nearest_hits <- distanceToNearest(
    peak_ranges,
    gene_ranges,
    ignore.strand = TRUE
  )

  # indices of peaks and their nearest genes
  peak_indices <- queryHits(nearest_hits)
  gene_indices <- subjectHits(nearest_hits)

  gene_names <- as.character(
    mcols(gene_ranges)$gene_name[gene_indices]
  )

  gene_ids <- as.character(
    mcols(gene_ranges)$gene_id[gene_indices]
  )

  # Use Ensembl ID when no gene symbol is available
  use_gene_id <- (
    is.na(gene_names) |
    gene_names == ""
  )

  gene_names[use_gene_id] <- gene_ids[use_gene_id]

  peaks$nearest_gene[peak_indices] <- gene_names

  peaks$distance_to_gene[peak_indices] <- as.integer(
    mcols(nearest_hits)$distance
  )
}

# ============================================================
# Write the peak annotation table
# ============================================================

# temp output file defs
output_peak_file <- file.path(workdir, "trait_ouptut_peaks.tsv")

peak_output_table <- peaks[
  ,
  .(
    trait = trait,
    snp,
    chromosome = chrom,
    position = pos,
    negLog10p = negLog10p,
    nearest_gene,
    distance_to_gene
  )
]

fwrite(
  peak_output_table,
  file = output_peak_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Peak table written to: ", output_peak_file)

# ============================================================
# Prepare CMplot input
# ============================================================

# replace with trait_df?

# Only send non-missing labels to CMplot
labeled_peaks <- peaks[
  !is.na(nearest_gene) &
  nearest_gene != ""
]

print('head of labeled_peaks:')
head(labeled_peaks)

cmplot_arguments <- list(
  Pmap = trait_dt,

  # Whole-genome Manhattan plot
  plot.type = "m",

  cex = 0.4,

  # Values are already -log10(P)
  LOG10 = FALSE,

  threshold = significance_threshold,
  threshold.col = "red",
  threshold.lty = 2,
  threshold.lwd = 1,
  amplify=FALSE,

  # Alternating chromosome colors
  col = c("#2166AC", "#67A9CF"),

  # Chromosomes 23, 24, and 25 are X, Y, and MT
  chr.labels =as.character(1:22),
  chr.labels.angle = 45,

  file.output = FALSE,
  verbose = FALSE
)

if (nrow(labeled_peaks) > 0L) {
  cmplot_arguments$highlight <- labeled_peaks$snp
  cmplot_arguments$highlight.col <- "red"
  cmplot_arguments$highlight.cex <- 0.8
  cmplot_arguments$highlight.text <- labeled_peaks$nearest_gene
  cmplot_arguments$highlight.text.cex <- 1
}



# ============================================================
# Create the plot
# ============================================================

# temp output file defs
output_plot_file <- file.path(workdir, "manat_plot.png")

# Check if the folder exists; if not, create it
#if (!dir.exists(dirname(output_plot_file))) {
#  dir.create(dirname(output_plot_file), recursive = TRUE)
#}

png(
  filename = output_plot_file,
  width = 4200,
  height = 2100,
  res = 300
)

# Increase the left margin (the second value in 'mar', default is usually 5.1)
#par(mar = c(bottom, left, top, right))  Default values: c(5.1, 4.1, 4.1, 2.1)
par(mar = c(4, 6, 2, 0))

plot_completed <- FALSE

tryCatch(
  {
    do.call(CMplot, cmplot_arguments)
    plot_completed <- TRUE
  },
  finally = {
    dev.off()
  }
)

if (!plot_completed || !file.exists(output_plot_file)) {
  stop("CMplot did not successfully create: ", output_plot_file)
}

message("Manhattan plot written to: ", output_plot_file)

# --- clean up workdir ---
#unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('DID NOT removed working directory:',workdir),quote=FALSE)
