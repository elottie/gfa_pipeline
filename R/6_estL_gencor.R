library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)
library(data.table)

# sample size affects genetic covariance and h2 but not intercept or genetic correlation
z_files = unlist(snakemake@input[["Z"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

# --- source helpful funcs ---
# make awks, sorts, and joins consistent across users
Sys.setenv(LC_ALL = "C")

harmon_helper_path <- "R/harmon_helpers.R"
ld_ref_helper_path <- "R/ld_ref_helpers.R"
# eventually need to switch to this
#helper_path <- "R/harmon_helpers.R"
source(harmon_helper_path)
source(ld_ref_helper_path)

# --- temp workdir for testing cleanliness -
workdir <- paste0("6_workdir_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = ""))
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

# temp output file defs
snps_in_ref_file <- file.path(workdir, "snps_in_ld_file.tsv")

# --- other random input setup ---
# if M is num of variants used to compute ld scores, should be constant?
M <- purrr:::map(1, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

# --- get ld ref across chromsomes ---
# previously would build ld ref across chr and use that to pull out just the matching snps from trait file
# but now, we are using snps that came out of gfa, so they are ld-pruned.  so ld-pruned snps << ld ref snps << trait snps
# so we could just leave the ld-pruned snps, but maybe save a little bit by doing intersection of ld-pruned and ref snps and pulling those out.  so I'll try
concat_ld_ref <- build_concat_ld_ref(ld_files = ld_files,
                                     snps_in_ref_file = snps_in_ref_file,
                                     delim="\t")

concat_status <- system(concat_ld_ref)
if (concat_status != 0) {
  stop("Failed to create unsorted reference with header")
}

print('completed ld ref concatenation')

# --- add the join between the ld-pruned and ld ref snps to reduce the number we pull out ---
ld_snps <- fread(snps_in_ref_file, header = TRUE, select = 1)[[1]]
head(ld_snps)
l2 <- as.numeric(scan(pipe(sprintf("awk 'NR > 1 {print $2}' %s", snps_in_ref_file)), what="character"))

# 1. Read headers only from all files
file_headers <- lapply(z_files, function(f) {
  names(fread(f, nrows = 0))
})

# 2. Check all files have exactly the same column names, in the same order
same_headers <- vapply(
  file_headers,
  function(hdr) identical(hdr, file_headers[[1]]),
  logical(1)
)

if (!all(same_headers)) {
  stop(
    "Not all z_files have the same column names. Problem files: ",
    paste(basename(z_files[!same_headers]), collapse = ", ")
  )
}

# 3. Extract z columns from file 1 only
z_cols <- file_headers[[1]][endsWith(file_headers[[1]], ".z")]

if (length(z_cols) == 0L) {
  stop("No columns ending in .z were found in the first file.")
}

# Preallocate output matrix
Z_hat <- matrix(
  NA_real_,
  nrow = length(ld_snps),
  ncol = length(z_cols),
  dimnames = list(ld_snps, str_replace(z_cols, ".z$", ""))
)

# Fill Z_hat matrix one file at a time

# make sure snps not matched across chr while processing files in loop
# don't think super necessary
prev_matched <- rep(FALSE, length(ld_snps))

for (j in seq_along(z_files)) {
  f <- z_files[[j]]

  z_file_filt <- fread(
    f,
    select = c("snp", z_cols)
  )

  # Direct row alignment. No duplicated SNP handling needed.
  ld_snp_index <- match(ld_snps, z_file_filt$snp)
  # keep previously match snps
  matched <- !is.na(ld_snp_index)

  snp_overlap_across_chr <- matched & prev_matched
  if (any(snp_overlap_across_chr)) {
    stop(
      "Some SNPs were matched in more than one z_file. Problem file: ",
      basename(f),
      ". Example SNPs: ",
      paste(head(ld_snps[snp_overlap_across_chr], 10), collapse = ", ")
    )
  }

  for (col_index in seq_along(z_cols)) {
    col <- z_cols[col_index]
    Z_hat[matched, col_index] <- z_file_filt[[col]][ld_snp_index[matched]]
  }

  prev_matched[matched] <- TRUE

  rm(z_file_filt, ld_snp_index)
  if (j %% 2 == 0) gc()
}

print(head(Z_hat))
print(tail(Z_hat))
print(dim(Z_hat))
print(sum(is.na(Z_hat)))

# --- obtain genetic correlation ---
R <- R_ldsc(Z_hat = Z_hat,
            ldscores = l2,
            ld_size = M,
            N = setNames(rep(1, ncol(Z_hat)), colnames(Z_hat)),
            return_gencov = TRUE,
            make_well_conditioned = FALSE # not needed we only need Rg
)

# --- save output ---
saveRDS(R, file=out)

# --- clean up workdir ---
unlink(workdir, recursive = TRUE, force = TRUE)
print(paste('removed working directory:',workdir),quote=FALSE)

