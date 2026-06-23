# --- functions for harmonizing (get correct Zs) on the fly ---

# extract only the columns we need from data files
map_gwas_info_cols <- function(gwas_info, trait, trait_col = "name") {
  
  row <- gwas_info[gwas_info[[trait_col]] == trait, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected 1 row for trait=", trait, " found ", nrow(row))

  as.list(row[1, c("snp","beta_hat","se","A1","A2","allele_freq","sample_size","chrom")])
}

# do the data harmonization
# gwas_info should be a data.table
# now returns rows in order of the snps_in_ref_file
# now returns all snps in ref file, with NAs for snps that are in ref file but not trait file
harmon_dat <- function(gwas_info, trait, snps_in_ref_file, return_ss=FALSE, return_alleles=FALSE, needs_invalid_snp_rm=FALSE) {
  
  full_trait <- gwas_info[name==trait, 'raw_data_path']
  print(paste("... processing trait:", full_trait))

  # make awks, sorts, and joins consistent across users
  Sys.setenv(LC_ALL = "C")
  # select rows we need.  sort/join would be more efficient but this is not huge for 1.25 mil line refs
  awk_trait_in_ref <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || ($1 in snps)' %s -",
    shQuote(full_trait), shQuote(snps_in_ref_file)
  )
  
  # read trait, just rows and columns we need
  dat_cols_sel <- map_gwas_info_cols(gwas_info, trait)
  trait_in_ref <- fread(cmd = awk_trait_in_ref, sep = "\t", header = TRUE,
            select = unname(unlist(dat_cols_sel)))  # pick what you need
  setnames(trait_in_ref, old = unname(unlist(dat_cols_sel)), new = names(dat_cols_sel))
  # grab pub ss, effect_or separately bc values, not column name
  pub_ss_val <- gwas_info[gwas_info[["name"]] == trait, "pub_sample_size", drop = FALSE]
  effect_or_flag <- gwas_info[gwas_info[["name"]] == trait, "effect_is_or", drop = FALSE]
  print(paste('pub_ss_val:',pub_ss_val))
  print(paste('effect_or_flag:',effect_or_flag))
  
  # now we want to make snps_in_ref the output table and join the trait file to it
  # so all snps_in_ref are in the final output, not all snps in trait file
  # read in ref file and deduplicate
  # header being true now means all files must have header even if just list of snps
  snps_in_ref <- fread(snps_in_ref_file, header = TRUE, select = 1)
  setnames(snps_in_ref, "snp")
  snps_in_ref <- snps_in_ref[, .N, by = "snp"][N == 1, .(snp)]
  # also deduplicate trait file so join does not overwrite rows
  print(paste('before removing dup snps before join with ref:',nrow(trait_in_ref)))
  dup_trait_snps <- trait_in_ref[duplicated(snp), unique(snp)]
  trait_in_ref <- trait_in_ref[!(snp %in% dup_trait_snps)]
  print(paste('after removing dup snps before join with ref:',nrow(trait_in_ref)))
  # select columns that will be added from trait file
  trait_cols <- setdiff(names(trait_in_ref), "snp")
  # in-place merge, mget gets names of columns added from trait file
  # so mem used is about trait_in_ref + merged_table
  snps_in_ref[
    trait_in_ref,
    on = "snp",
   (trait_cols) := mget(paste0("i.", trait_cols))
  ]
  # does not create new object, just new name for same object
  filt_trait <- snps_in_ref
  # cleanup (removal of snps_in_ref not necessary, just removing name)
  rm(dup_trait_snps, trait_in_ref, snps_in_ref)
  gc()

  print('head filt_trait after reading in:')
  print(head(filt_trait))

  if (toupper(effect_or_flag) %in% c("TRUE", "T")) {
    filt_trait[, beta_hat := {
      filt_trait[, beta_hat := as.numeric(beta_hat)]
      # doing in 2 steps like this slightly more optimal for avoiding warnings fifelse would produce from doing log first
      filt_trait[is.na(beta_hat) | beta_hat <= 0, beta_hat := NA_real_]
      filt_trait[!is.na(beta_hat), beta_hat := log(beta_hat)]
    }]
  }

  print('head filt_trait after standardizing betas to log(OR):')
  print(head(filt_trait))

  if (needs_invalid_snp_rm){
    gwas_format(filt_trait,
		snp="snp", beta_hat="beta_hat", se="se", A1="A1", A2="A2",
                chrom="chrom",sample_size="sample_size", allele_freq="allele_freq",compute_pval=FALSE)
  } else {
    GFA:::align_beta(filt_trait)
  }
  print('head filt_trait_harmon:')
  print(head(filt_trait))

  # we ensured beta_hat was numeric during align_beta.  just guard se
  filt_trait[, Z := {
    se_num <- as.numeric(se)
    fifelse(!is.na(se_num) & se_num != 0, beta_hat / se_num, NA_real_)
  }]
  print('head filt_trait with z:')
  print(head(filt_trait))

  # fill in any NAs in sample size w/ pub ss
  print(nrow(filt_trait[is.na(sample_size) | sample_size == "",]))
  filt_trait[is.na(sample_size) | sample_size == "", sample_size := pub_ss_val]
  print('head filt_trait after filling in missing ss:')
  print(head(filt_trait))

  if (return_ss) {
    return(list(snps = filt_trait[["snp"]],
		Z = filt_trait[["Z"]],
                ss = filt_trait[["sample_size"]]))

  } else if (return_alleles) {
    return(list(snps = filt_trait[["snp"]],
		Z = filt_trait[["Z"]],
		ref = filt_trait[["A2"]],
		alt = filt_trait[["A1"]]))
  } else {
    return(list(snps = filt_trait[["snp"]],
                Z = filt_trait[["Z"]]))
  }
}
