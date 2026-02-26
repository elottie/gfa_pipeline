library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)


# sample size affects genetic covariance and h2 but not intercept or genetic correlation
gwas_info <- read_csv(snakemake@input[["gwas_info"]])
strip_list <- readRDS(snakemake@input[["strip_file"]])
strip_number <- as.numeric(snakemake@wildcard[["strip_number"]])
out <- snakemake@output[["out"]]
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])

## make blocks data frame
nblocks <- length(strip_list)

results <- list()
for(s2 in strip_number:nblocks){
	if(s2 == strip_number){
		# read data for set 1
		z_files1 = strip_number
	}
	if(s2 ==strip_number + 1){
		# read set 2
	}else if(s2 > strip_number + 1){
		# drop old set 2, read new set 2
	}
	# do ldsc caclulations, store in results
}


saveRDS(ret, file=out)



