
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
  trait_memory_gb = 5,
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

  max_traits_per_block <- floor(memory_limit_gb / (blocks_at_once * trait_memory_gb))

  if (max_traits_per_block < min_traits_per_block) {
    stop(sprintf(
      "max_traits_per_block=%d < min_traits_per_block=%d given memory_limit_gb=%s, blocks_at_once=%s, trait_memory_gb=%s.",
      max_traits_per_block, min_traits_per_block, memory_limit_gb, blocks_at_once, trait_memory_gb
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
