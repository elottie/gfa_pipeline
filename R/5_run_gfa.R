library(dplyr)
library(purrr)
library(stringr)
library(GFA)


out <- snakemake@output[["out"]]
Z_file <- snakemake@input[["Z"]]
R_est_file <- snakemake@input[["R"]]
params_file <- snakemake@params[["params_file"]]
max_snp <- as.numeric(snakemake@params[["max_snps"]])
seed <- snakemake@wildcards[["fs"]]

load(Z_file)  # RData from make_nice_data
print(R_est_file)
R <- readRDS(R_est_file)

set.seed(seed)


if(params_file == "default"){
  params <- gfa_default_parameters()
}else{
  params <- readRDS(params_file)
}


# move this to make_nice_data.R? ---
# I think keep here and the slightly inefficient z_order.  keeping here avoids passing ton of params to make_nice_data
# z_order seemed not to be a big problem at least w/ 100s of traits.  may rewrite if becomes an issue

# Read in data
# we are not doing max_snps for ss here!!
if(nrow(Z_hat) > max_snp){
    ix <- sample(seq(nrow(Z_hat)), size = max_snp, replace = FALSE)
    Z_hat <- Z_hat[ix,]
}

snps <- rownames(Z_hat)
nms_z <- colnames(Z_hat)

# if(str_ends(R_est_file, "none_R.txt")){
#   R <- list(names = nms, R = diag(length(nms), nrow = ntrait))
# }else{
stopifnot(identical(rownames(R),colnames(R)))
nms_r <- colnames(R)
stopifnot(all(nms_r %in% nms_z))
z_order <- match(nms_r, nms_z)  # matches order of nms_r
ss <- ss[,z_order]
Z_hat <- Z_hat[,z_order]
#R$R <- cov2cor(R$R)
#}

rownames(R) <- colnames(R) <- NULL


#N <- apply(SS, 2, median)
# ---

t <- system.time(f <- gfa_fit(Z_hat = Z_hat,
                                N = ss,
                                R = R,
                                params = params 
                                #mode = "z-score",
                                #method = "fixed_factors"
				)
		)


f$snps <- snps
f$names <- nms_r
f$time <- t
saveRDS(f, file=out)

