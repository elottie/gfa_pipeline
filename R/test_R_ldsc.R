library(GFA)

# --- make toy setup to test Jean's updated R_ldsc ---

# Toy Z-score matrix: rows = SNPs, columns = traits
Z_hat <- matrix(c(
  1.2, 0.8, 0.9, 1.9,  # SNP1
  1.1, 0.7, 1.5, 0.2,  # SNP2
  0.9, 0.5, 1.8, 0.7,  # SNP3
  1.0, 0.6, 2.1, 0.9    # SNP4
), nrow = 4, byrow = TRUE)

# LD scores (numeric vector) for each SNP.  come from LD ref panel, not me
ldscores <- c(1.5, 2.0, 1.8, 2.2)

# LD size, num of variants used to compute LD scores, right now arbitrary scalar
ld_size <- 4

#N: sample sizes per trait (if equal across SNPs, can be a vector)
N <- c(1000, 1200, 1200, 1100)
blocks <- NULL
ncores <- 1

# the new option, comparisons of traits.  comparisons MUST be df, not matrix
comparisons <- data.frame(trait1 = c(1,3),
                          trait2 = c(2,4))
print(comparisons)

result <- R_ldsc(
  Z_hat = Z_hat,
  ldscores = ldscores,
  ld_size = ld_size,
  N = N,
  blocks = blocks,
  ncores = ncores,
  make_well_conditioned = FALSE,
  comparisons = comparisons
)
print(result)
