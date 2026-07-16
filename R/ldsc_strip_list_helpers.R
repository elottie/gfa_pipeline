
assign_traits <- function(traits, nsets) {
  n <- length(traits)
  base_num_traits <- n %/% nsets
  bigger_sets <- n %% nsets

  out <- vector("list", nsets)
  idx <- 1

  for (i in seq_len(nsets)) {
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
  max_traits_per_set = 20,
  sets_at_once = 2,
  min_traits_per_set = 2
) {
  if (!(name_col %in% names(gwas_info))) {
    stop(sprintf("Column '%s' not found in %s.", name_col, gwas_info_path))
  }

  traits <- as.character(gwas_info[[name_col]])
  traits <- traits[!is.na(traits)]
  ntraits <- length(traits)

  if (ntraits < min_traits_per_set) {
    stop(sprintf("Need at least %d traits; got %d.", min_traits_per_set, ntraits))
  }

  # got this from binary search & lin reg, hardcoded.  this is slope for each additional total trait analyzed at once
  # add safety factor of 2 
#  safety_fac <- 5
#  traits_slope <- 40*safety_fac / 1024

#  max_traits_per_block <- floor(mem_for_traits / traits_slope)
#  print(paste('max traits per block from mem avail for traits:',max_traits_per_block),quote=F)
  print(paste('max traits per set from user input (>40 results in runtimes >13h):',max_traits_per_set),quote=F)

  if (max_traits_per_set < min_traits_per_set) {
    stop(sprintf(
      "max_traits_per_set=%d < min_traits_per_set=%d given sets_at_once=%s.",
      max_traits_per_set, min_traits_per_set, sets_at_once
    ))
  }

  # Feasible nblocks range
  nsets_min <- ceiling(ntraits / max_traits_per_set)  # enough blocks to keep size <= max
  nsets_max <- ntraits %/% min_traits_per_set         # not too many blocks so size >= min

  if (nsets_min > nsets_max) {
    stop(sprintf(
      "Impossible to satisfy block size constraints: ntraits=%d, min=%d, max=%d.",
      ntraits, min_traits_per_set, max_traits_per_set
    ))
  }

  # Fewest blocks (goal)
  nsets <- nsets_min

  strip_list <- assign_traits(traits, nsets)

  # Validate sizes
  sizes <- vapply(strip_list, length, integer(1))
  if (min(sizes) < min_traits_per_set || max(sizes) > max_traits_per_set) {
    stop(sprintf(
      "Internal error: set sizes [%s] violate [%d, %d].",
      paste(sizes, collapse = ", "),
      min_traits_per_set, max_traits_per_set
    ))
  }

  strip_list
}

# Example:
# blocks <- make_trait_sets("traits.csv", name_col = "name")
# str(blocks)
