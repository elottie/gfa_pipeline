library(readr)
library(dplyr)

comparison_files <- list.files("ldsc_block", pattern = "set_\\d+_\\d+\\.csv", full.names = TRUE)
print(comparison_files)

# Helper to load trait IDs for a set (from comparison file)
get_trait_ids <- function(comparison_file, setnum) {
  comparisons <- read_csv(comparison_file)
  setnums <- as.numeric(regmatches(basename(comparison_file), gregexpr("\\d+", basename(comparison_file)))[[1]])
  if(setnum == setnums[1]) {
    unique(comparisons$trait1)
  } else {
    unique(comparisons$trait2)
  }
}

make_trait_df <- function(trait_ids, set_name) {
  data.frame(trait = trait_ids, value = runif(length(trait_ids)), set = set_name)
}

set_dfs <- list() # to hold only what's needed

comparison_results <- data.frame(
  s1 = character(),
  s2 = character(),
  trait1 = character(),
  trait2 = character(),
  se = numeric(),      # residual covariance
  sg = numeric(),      # genetic covariance
  rg = numeric(),      # genetic correlation
  stringsAsFactors = FALSE
)

for (i in seq_along(comparison_files)) {
  file <- comparison_files[i]
  setnums <- as.numeric(regmatches(basename(file), gregexpr("\\d+", basename(file)))[[1]])
  print(setnums)
  s1name <- paste0("s", setnums[1])
  s2name <- paste0("s", setnums[2])

  print(paste("s1name",s1name))
  print(paste("s2name",s2name))

  # Only create set dataframe if it doesn't already exist/nor just deleted
  if (!s1name %in% names(set_dfs)) {
    trait_ids_s1 <- get_trait_ids(file, setnums[1])
    set_dfs[[s1name]] <- make_trait_df(trait_ids_s1, s1name)
  }
  if (!s2name %in% names(set_dfs)) {
    trait_ids_s2 <- get_trait_ids(file, setnums[2])
    set_dfs[[s2name]] <- make_trait_df(trait_ids_s2, s2name)
  }

  comparisons <- read_csv(file)
  # head(comparisons) # Remove unless you want to print each time

  cat("\n--- Comparing traits in", s1name, "to traits in", s2name, "---\n")

  for (j in 1:nrow(comparisons)) {
    t1 <- comparisons$trait1[j]
    t2 <- comparisons$trait2[j]
    cat(sprintf("Comparing %s (from %s) to %s (from %s)\n", t1, s1name, t2, s2name))
    # HERE would be equivalent of ldsc_rg.  potentially gecko version
    #result_vec <- ldsc_rg(...)  
    # not sure these are the actual outputs
    #result_vec <- c(
    #    int = runif(1, 0, 0.1),
    #    int_se = runif(1, 0, 0.01),
    #    h2 = runif(1, 0, 1),
    #    h2_se = runif(1, 0, 0.1)
    #)

    result_vec <- c(
      se = runif(1, 0, 1),           # Dummy residual covariance
      sg = runif(1, 0, 1),           # Dummy genetic covariance
      rg = runif(1, -1, 1)           # Dummy genetic correlation (could range -1 to 1)
    )
    comparison_results <- rbind(comparison_results,
        data.frame(
            s1 = s1name,
            s2 = s2name,
            trait1 = t1,
            trait2 = t2,
            se = result_vec[1],
            sg = result_vec[2],
            rg = result_vec[3],
            stringsAsFactors = FALSE
        )
    )

  }

  # ---- TRIANGLE PROCESSING ----
  # For set_1_2: also process t1, t2
  # For set_1_3: also process t3
  # For set_1_4: also process t4
  triangle_dir <- "ldsc_block" # change if triangles are elsewhere

  if (setnums[2] == 2) { # For set_1_2, do t1 AND t2
    print('in setnums2')
    for (tri_num in 1:2) {
      tri_file <- file.path(triangle_dir, paste0("triangle_", tri_num, ".csv"))
      print(paste('tri_file:',tri_file))
      if (file.exists(tri_file)) {
        print('tri_file exists')
        triangle_comparisons <- read_csv(tri_file)
        cat(sprintf("\nProcessing triangle comparisons in %s:\n", basename(tri_file)))
        for (k in 1:nrow(triangle_comparisons)) {
          tri_t1 <- triangle_comparisons$trait1[k]
          tri_t2 <- triangle_comparisons$trait2[k]
          cat(sprintf("Triangle: Comparing %s to %s\n", tri_t1, tri_t2))
        }
      }
    }
  } else if (setnums[2] >= 3) { # For set_1_3, set_1_4, etc: just t3, t4, ...
    tri_file <- file.path(triangle_dir, paste0("triangle_", setnums[2], ".csv"))
    if (file.exists(tri_file)) {
      triangle_comparisons <- read_csv(tri_file)
      cat(sprintf("\nProcessing triangle comparisons in %s:\n", basename(tri_file)))
      for (k in 1:nrow(triangle_comparisons)) {
        tri_t1 <- triangle_comparisons$trait1[k]
        tri_t2 <- triangle_comparisons$trait2[k]
        cat(sprintf("Triangle: Comparing %s to %s\n", tri_t1, tri_t2))
      }
    }
  }

  # ---- END TRIANGLE PROCESSING ----

  # Drop s2 after comparison if not the last set (you could adapt this logic as needed)
  set_dfs[[s1name]] <- NULL # optional, if you only want to keep s2 (or s3...) for next round
  cat(sprintf("Dropped dataframe for %s\n", s2name))
}

cat("\n--- All comparisons (and triangles) complete! ---\n")

cat('\n--- combining results (questionable) ---\n')

library(tidyr)

# Se_mat (Residual covariance)
Se_mat <- comparison_results %>%
    select(trait1, trait2, se) %>%
    pivot_wider(names_from = trait2, values_from = se)

rownames(Se_mat) <- Se_mat$trait1
Se_mat <- Se_mat[, !(names(Se_mat) == "trait1")]
Se_mat <- as.matrix(Se_mat)

# Sg_mat (Genetic covariance)
Sg_mat <- comparison_results %>%
    select(trait1, trait2, sg) %>%
    pivot_wider(names_from = trait2, values_from = sg)

rownames(Sg_mat) <- Sg_mat$trait1
Sg_mat <- Sg_mat[, !(names(Sg_mat) == "trait1")]
Sg_mat <- as.matrix(Sg_mat)

# Rg_mat (Genetic correlation)
Rg_mat <- comparison_results %>%
    select(trait1, trait2, rg) %>%
    pivot_wider(names_from = trait2, values_from = rg)

rownames(Rg_mat) <- Rg_mat$trait1
Rg_mat <- Rg_mat[, !(names(Rg_mat) == "trait1")]
Rg_mat <- as.matrix(Rg_mat)

cat('\n--- all results combined! ---\n')

cat('\n--- storing results in a list ---\n')
output <- list(
    Se = Se_mat,
    Sg = Sg_mat,
    Rg = Rg_mat
    # add Ve, Vg, VRg
)
cat('\n--- results stored in a list like so: ---\n')
str(output)
