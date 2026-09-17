suppressPackageStartupMessages({
  library(data.table)
  library(CMplot)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

plot_manhattan <- function(dt,
			   analysis_name,
			   genome_build,
			   sig_thresh_negLog10p,
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

  sig_thresh_negLog10p <- as.numeric(sig_thresh_negLog10p)
  peak_window <- as.numeric(peak_window)

  message("Analysis name: ", analysis_name)
  print("Head of dt to plot:")
  print(head(dt))
  message("Genome build: ", genome_build)
  message("Significance threshold, -log10(p): ", sig_thresh_negLog10p)
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
  
  # Keep a copy of every SNP above the significance threshold for output table, only peaks are in manhat plot
  significant_snps <- copy(
    dt[negLog10p >= sig_thresh_negLog10p]
  )

  # Temporary row ID allows assignments to survive candidate reordering.
  significant_snps[, significant_row_id := .I]
  significant_snps[, peak_id := NA_integer_]

  # peak_candidates will be reduced during peak selection, while
  # significant_snps remains complete.
  peak_candidates <- copy(significant_snps)
  
  # Start with the strongest association
  setorder(peak_candidates, -negLog10p)

  peaks <- peak_candidates[0]

  while (nrow(peak_candidates) > 0L) {

    lead_snp <- copy(peak_candidates[1])
    peak_number <- nrow(peaks) + 1L

    # Identify all currently unassigned significant SNPs belonging to
    # this peak. This follows the same rule used to select the peaks.
    peak_member_snps <- peak_candidates[
      chrom == lead_snp$chrom &
        abs(pos - lead_snp$pos) <= peak_window,
      significant_row_id
    ]

    # Record the peak associated with each significant SNP.
    significant_snps[
      peak_member_snps,
      peak_id := peak_number
    ]

    # Give the lead SNP the same peak identifier.
    lead_snp[, peak_id := peak_number]

    peaks <- rbind(
      peaks,
      lead_snp,
      use.names = TRUE
    )

    # Remove significant peak_candidates assigned to this peak.
    peak_candidates <- peak_candidates[
      chrom != lead_snp$chrom |
        abs(pos - lead_snp$pos) > peak_window
    ]
  }

  # Renumber peaks according to their genomic order in the original dt.
  peak_id_map <- peaks[
    order(significant_row_id),
    .(old_peak_id = peak_id)
  ]
  peak_id_map[, new_peak_id := .I]

  # Apply the new IDs to all significant SNPs.
  significant_snps[
    peak_id_map,
    on = .(peak_id = old_peak_id),
    peak_id := i.new_peak_id
  ]

  # Apply the same IDs to the lead-peak table.
  peaks[
    peak_id_map,
    on = .(peak_id = old_peak_id),
    peak_id := i.new_peak_id
  ]

  setorder(peaks, peak_id)

  message("Significant peaks selected for labeling: ", nrow(peaks))
  message("Total significant SNPs for table: ", nrow(significant_snps))

  # ============================================================
  # Find the nearest gene for every significant SNP
  # ============================================================

  significant_snps[, `:=`(
    nearest_gene = NA_character_,
    distance_to_gene = NA_integer_
  )]

  if (nrow(significant_snps) > 0L) {
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

    # Create genomic ranges for every significant SNP.
    significant_ranges <- GRanges(
      seqnames = significant_snps$chrom,
      ranges = IRanges(
        start = as.integer(significant_snps$pos),
        end = as.integer(significant_snps$pos)
      )
    )
    gene_ranges <- genes(
      annotation_database,
      columns = c("gene_id", "gene_name")
    )

    nearest_hits <- distanceToNearest(
      significant_ranges,
      gene_ranges,
      ignore.strand = TRUE
    )

    # Indices of significant SNPs and their nearest genes.
    snp_indices <- queryHits(nearest_hits)
    gene_indices <- subjectHits(nearest_hits)
    gene_names <- as.character(
      mcols(gene_ranges)$gene_name[gene_indices]
    )
    gene_ids <- as.character(
      mcols(gene_ranges)$gene_id[gene_indices]
    )
    # Use the Ensembl ID when no gene symbol is available.
    use_gene_id <- is.na(gene_names) | gene_names == ""
    gene_names[use_gene_id] <- gene_ids[use_gene_id]

    significant_snps[
      snp_indices,
      nearest_gene := gene_names
    ]

    significant_snps[
      snp_indices,
      distance_to_gene := as.integer(mcols(nearest_hits)$distance)
    ]
  }

  # ============================================================
  # Recover the annotated lead SNPs for Manhattan-plot labeling
  # ============================================================

  # match() preserves the original peak-selection order.
  peaks <- significant_snps[
    match(
      peaks$significant_row_id,
      significant_snps$significant_row_id
    )
  ]

  # ============================================================
  # Write every significant SNP to the output table
  # ============================================================

  significant_output <- significant_snps[
    ,
    .(
      analysis = analysis_name,
      snp,
      chromosome = chrom,
      position = pos,
      negLog10p,
      peak_id,
      nearest_gene,
      distance_to_gene
    )
  ]

  fwrite(
    significant_output,
    file = peaks_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  message("All significant SNPs and nearest genes written to: ",peaks_file)  

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
    threshold = sig_thresh_negLog10p,
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
