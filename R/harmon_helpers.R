# --- functions for harmonizing (get correct Zs) on the fly ---

# extract only the columns we need from data files
map_gwas_info_cols <- function(gwas_info, trait, trait_col = "name") {
  
  row <- gwas_info[gwas_info[[trait_col]] == trait, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected 1 row for trait=", trait, " found ", nrow(row))

  as.list(row[1, c("snp","beta_hat","se","A1","A2","af","sample_size")])
}

# do the data harmonization
# gwas_info should be a data.table
harmon_dat <- function(gwas_info, trait, snps_in_ref_file, return_ss=FALSE) {
  
  full_trait <- gwas_info[name==trait, 'raw_data_path']

  print(paste("... processing trait:", full_trait))

  awk_trait_in_ref <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || ($1 in snps)' %s -",
    shQuote(full_trait), shQuote(snps_in_ref_file)
  )

  dat_cols_sel <- map_gwas_info_cols(gwas_info, trait)
  filt_trait <- fread(cmd = awk_trait_in_ref, sep = "\t", header = TRUE,
            select = unname(unlist(dat_cols_sel)))  # pick what you need
  setnames(filt_trait, old = unname(unlist(dat_cols_sel)), new = names(dat_cols_sel))
  # grab pub ss separately bc it's just a value, not column name
  pub_ss_val <- gwas_info[gwas_info[["name"]] == trait, "pub_sample_size", drop = FALSE]
  effect_or_flag <- gwas_info[gwas_info[["name"]] == trait, "effect_is_or", drop = FALSE]
  print(paste('pub_ss_val:',pub_ss_val))
  print(paste('effect_or_flag:',effect_or_flag))
  # reorder rows to match snps_in_ref_file
  snps_in_ref <- fread(snps_in_ref_file, header = FALSE, select = 1)
  setnames(snps_in_ref, 1, "snp")
  snps_in_ref <- unique(snps_in_ref, by = "snp")
  filt_trait <- filt_trait[snps_in_ref, on = "snp", nomatch = 0]

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

  # fix the readability of making beta_hat numeric
  GFA:::align_beta(filt_trait)
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

  if(return_ss){
    return(list(Z = filt_trait[["Z"]],
                ss = filt_trait[["sample_size"]]))

  } else{
    return(list(Z = filt_trait[["Z"]]))
  }
}
