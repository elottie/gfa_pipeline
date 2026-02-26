get_set <- function(mat, block_size, i, j, trait_names) {
  row_start <- (i - 1) * block_size + 1
  row_end   <- min(i * block_size, nrow(mat))
  col_start <- (j - 1) * block_size + 1
  col_end   <- min(j * block_size, ncol(mat))
  block <- mat[row_start:row_end, col_start:col_end]
  row_names <- trait_names[row_start:row_end]
  col_names <- trait_names[col_start:col_end]
  list(block=block, row_names=row_names, col_names=col_names)
}

get_triangle <- function(mat, block_size, i, trait_names) {
  set <- get_set(mat, block_size, i, i, trait_names)
  mask <- upper.tri(set$block)
  inds <- which(mask, arr.ind=TRUE)
  data.frame(
    trait1 = set$row_names[inds[,1]],
    trait2 = set$col_names[inds[,2]],
    value  = set$block[mask]
  )
}

process_strips_write <- function(mat, block_size = 1000, out_dir = "results") {
  n <- nrow(mat)
  stopifnot(n == ncol(mat))
  trait_names <- rownames(mat)
  if (is.null(trait_names)) trait_names <- as.character(1:n)
  num_blocks <- ceiling(n / block_size)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (i in 1:num_blocks) {
    # Indices for this block
    row_start <- (i - 1) * block_size + 1
    row_end   <- min(i * block_size, n)
    row_names <- trait_names[row_start:row_end]

    print(paste("rowstart:",row_start))
    print(paste("rowend:",row_end))
    
    # Diagonal block: within-block upper triangle
    block <- mat[row_start:row_end, row_start:row_end]
    inds <- which(upper.tri(block), arr.ind=TRUE)
    triangle_df <- data.frame(
      trait1 = row_names[inds[,1]],
      trait2 = row_names[inds[,2]],
      value  = block[upper.tri(block)]
    )
    # Save triangle to disk
    triangle_file <- file.path(out_dir, sprintf("triangle_%d.csv", i))
    write.csv(triangle_df, triangle_file, row.names = FALSE)
    
    # Process all blocks above the diagonal (off-diagonal blocks)
    for (j in (i+1):num_blocks) {
      if (j < row_start | j > num_blocks) next # don't process invalid blocks
      col_start <- (j - 1) * block_size + 1
      col_end   <- min(j * block_size, n)
      col_names <- trait_names[col_start:col_end]
      set_block <- mat[row_start:row_end, col_start:col_end]

      print(paste("colstart:",col_start))
      print(paste("colend:",col_end))
      
      inds_set <- expand.grid(row = 1:nrow(set_block), col = 1:ncol(set_block))
      set_df <- data.frame(
        trait1 = row_names[inds_set$row],
        trait2 = col_names[inds_set$col],
        value  = as.vector(set_block)
      )
      set_file <- file.path(out_dir, sprintf("set_%d_%d.csv", i, j))
      write.csv(set_df, set_file, row.names = FALSE)
      
      rm(set_block, set_df, inds_set)
      gc()
    }
    rm(block, triangle_df, inds)
    gc()
  }
}

#How to Run
traits <- read.csv("../5e5Sig_Herit_Mets_8ForLDSCStrip.csv", stringsAsFactors = FALSE)$name
# Now traits is a vector of trait names

#2. Create the Square Matrix
#A. Fill with NA or zeros
n <- length(traits)
mat <- matrix(NA, nrow = n, ncol = n,
              dimnames = list(traits, traits))

process_strips_write(mat, block_size = 2, out_dir = "./ldsc_block/")
