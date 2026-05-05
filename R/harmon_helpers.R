# --- functions for harmonizing (get correct Zs) on the fly ---

# extract only the columns we need from data files
map_gwas_info_cols <- function(gwas_info, trait, trait_col = "name") {
  
  row <- gwas_info[gwas_info[[trait_col]] == trait, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected 1 row for trait=", trait, " found ", nrow(row))

  as.list(row[1, c("beta_hat","se","A1","A2","af","sample_size")])
}

# do the data harmonization
# gwas_info should be a data.table
harmon_dat <- function(gwas_info, trait, snps_in_ref, return_ss=FALSE) {
  
  full_trait <- gwas_info[name==trait, 'raw_data_path']

  print(paste("... processing trait:", full_trait))

  awk_trait_in_ref <- sprintf(
    "zcat %s | awk 'NR==FNR {snps[$1]=1; next} FNR==1 || snps[$1]' %s -",
    shQuote(full_trait), shQuote(snps_in_ref)
  )

  dat_cols_sel <- map_gwas_info_cols(gwas_info, trait)
  filt_trait <- fread(cmd = awk_trait_in_ref, sep = "\t", header = TRUE,
            select = unname(unlist(dat_cols_sel)))  # pick what you need
  setnames(filt_trait, old = unname(unlist(dat_cols_sel)), new = names(dat_cols_sel))

  print('head filt_trait after reading in:')
  print(head(filt_trait))

#  filtered_trait_harmon <- GFA:::aligss_beta(filtered_trait,beta_name,af_name)
  GFA:::align_beta(filt_trait)
  print('head filt_trait_harmon:')
  print(head(filt_trait))

  # we ensured beta_hat was numeric during align_beta.  just guard se
  filt_trait[, se_num := as.numeric(se)]
  filt_trait[, Z := fifelse(se_num != 0, beta_hat / se_num, NA_real_)]
  print('head filt_trait with z:')
  print(head(filt_trait))
  
  if(return_ss){
    return(list(Z = filt_trait[["Z"]],
                ss = filt_trait[["sample_size"]]))

  } else{
    return(list(Z = filt_trait[["Z"]]))
  }
}
