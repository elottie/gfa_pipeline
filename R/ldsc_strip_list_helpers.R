
assign_traits <- function(traits, nblocks) {
  n <- length(traits)
  base_num_traits <- n %/% nblocks
  bigger_sets <- n %% nblocks

  out <- vector("list", nblocks)
  idx <- 1

  for (i in seq_len(nblocks)) {
    set_size <- base_num_traits + if (i <= bigger_sets) 1L else 0L
    if (set_size == 0L) {
      out[[i]] <- character(0)
    } else {
      out[[i]] <- traits[idx:(idx + set_size - 1L)]
    }
    idx <- idx + set_size
  }
  out
}

make_trait_sets <- function(
  gwas_info,
  name_col = "name",
  memory_limit_gb = 150,
  blocks_at_once = 2,
  min_traits_per_block = 2
) {
  if (!(name_col %in% names(gwas_info))) {
    stop(sprintf("Column '%s' not found in %s.", name_col, gwas_info_path))
  }

  traits <- as.character(gwas_info[[name_col]])
  traits <- traits[!is.na(traits)]
  ntraits <- length(traits)

  if (ntraits < min_traits_per_block) {
    stop(sprintf("Need at least %d traits; got %d.", min_traits_per_block, ntraits))
  }

  # R needs MB ish just to exist
  R_intercept <- 1200 / 1024

  mem_for_traits <- memory_limit_gb - R_intercept
  print(paste('mem avail for traits in gb:',mem_for_traits),quote=FALSE)

  # say ref panel is 1.25 million lines.  this can / should be checked.  anyway, by awks, this is the max # of rows for traits
  #L <- 1.25e6
  # we are storing 2 vectors for each trait:  the Z and the sample_size
  # so the size for one block:  bytes_1_block = 8 bytes * (2 vectors per trait * L rows * T traits in block) = 16 L T_in_block
  # the size for 2 blocks:  bytes_2_blocks = 2 * bytes_1_block = 32 L T_in_block
  # so max # of traits in a block:  T_in_block = bytes / (32 * L)

  # got this from binary search & lin reg, hardcoded.  this is slope for each additional total trait analyzed at once
  # add safety factor of 2 
  safety_fac <- 2
  traits_slope <- 35*safety_fac / 1024

  max_traits_per_block <- floor(mem_for_traits / traits_slope)
  print(paste('max traits per block from mem avail for traits:',max_traits_per_block),quote=F)

  if (max_traits_per_block < min_traits_per_block) {
    stop(sprintf(
      "max_traits_per_block=%d < min_traits_per_block=%d given memory_limit_gb=%s, blocks_at_once=%s.",
      max_traits_per_block, min_traits_per_block, memory_limit_gb, blocks_at_once
    ))
  }

  # Feasible nblocks range
  nblocks_min <- ceiling(ntraits / max_traits_per_block)  # enough blocks to keep size <= max
  nblocks_max <- ntraits %/% min_traits_per_block         # not too many blocks so size >= min

  if (nblocks_min > nblocks_max) {
    stop(sprintf(
      "Impossible to satisfy block size constraints: ntraits=%d, min=%d, max=%d.",
      ntraits, min_traits_per_block, max_traits_per_block
    ))
  }

  # Fewest blocks (goal)
  nblocks <- nblocks_min

  strip_list <- assign_traits(traits, nblocks)

  # Validate sizes
  sizes <- vapply(strip_list, length, integer(1))
  if (min(sizes) < min_traits_per_block || max(sizes) > max_traits_per_block) {
    stop(sprintf(
      "Internal error: block sizes [%s] violate [%d, %d].",
      paste(sizes, collapse = ", "),
      min_traits_per_block, max_traits_per_block
    ))
  }

  strip_list
}

# Example:
# blocks <- make_trait_sets("traits.csv", name_col = "name")
# str(blocks)
