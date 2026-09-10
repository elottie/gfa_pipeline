suppressPackageStartupMessages({
  library(data.table)
  library(CMplot)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

plot_manhattan <- function(dt,
			   analysis_name,
			   genome_build,
			   significance_threshold,
			   peak_window,
                           plot_file,
                           peaks_file) {

  # ============================================================
  # inputs, outputs, and parameters
  # ============================================================

  # expects dt has columns snp, chrom, pos, negLog10p.  we will find nearest genes for you	
  stopifnot(all(c('snp','chrom','pos','negLog10p') %in% colnames(dt)))

  # Check if the analysis folder exists; if not, create it
  if (!dir.exists(dirname(plot_file))) {
    dir.create(dirname(plot_file), recursive = TRUE)
  }
  if (!dir.exists(dirname(peaks_file))) {
    dir.create(dirname(peaks_file), recursive = TRUE)
  }

  significance_threshold <- as.numeric(significance_threshold)
  peak_window <- as.numeric(peak_window)

  message("Analysis name: ", analysis_name)
  print("Head of dt to plot:")
  print(head(dt))
  message("Genome build: ", genome_build)
  message("Significance threshold, -log10(p): ", significance_threshold)
  message("Peak window: ", peak_window, " bp")

  # ============================================================
  # Clean
  # ============================================================

  dt <- dt[chrom %between% c(1, 22)]

  # valid snp selection assumed handled earlier

  # CMplot uses variant identifiers for highlighting. Make them unique
  # rm all snps that ever occur more than once
  dt <- dt[, if (.N == 1L) .SD, by = snp]

  setorder(dt, chrom, pos)

  message("Valid variants: ", nrow(dt))
  message(
    "Chromosomes represented: ",
    paste(unique(dt$chrom), collapse = ", ")
  )

  # ============================================================
  # Select significant peaks for gene labeling
  # ============================================================
  #
  # Every valid variant is included in the Manhattan plot.
  # peak_window only determines which significant variants receive labels.

  candidates <- copy(
    dt[negLog10p >= significance_threshold]
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

  peak_output <- peaks[
    ,
    .(
      analysis = analysis_name,
      snp,
      chromosome = chrom,
      position = pos,
      negLog10p = negLog10p,
      nearest_gene,
      distance_to_gene
    )
  ]

  fwrite(
    peak_output,
    file = peaks_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  message("Manhattan peak table written to: ", peaks_file)

  # ============================================================
  # Prepare CMplot input
  # ============================================================

  # Only send non-missing labels to CMplot
  labeled_peaks <- peaks[
    !is.na(nearest_gene) &
    nearest_gene != ""
  ]

  cmplot_arguments <- list(
    Pmap = dt,

    # Whole-genome Manhattan plot
    plot.type = "m",

    # Alternating chromosome colors
    col = c("#2166AC", "#67A9CF"),

    # dot size
    cex = 0.4,

    # Chromosomes 23, 24, and 25 are X, Y, and MT
    chr.labels = as.character(1:22),
    chr.labels.angle = 45,
     
    # Values are already -log10(P)
    LOG10 = FALSE,

    # significance line
    threshold = significance_threshold,
    threshold.col = "red",
    threshold.lty = 2,
    threshold.lwd = 1,
    amplify=FALSE,

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

  png(
    filename = plot_file,
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

  if (!plot_completed || !file.exists(plot_file)) {
    stop("CMplot did not successfully create: ", plot_file)
  }

  message("Manhattan plot written to: ", plot_file)

}
