# add_zscore.awk
# expects -v beta_name="..." -v se_name="..."

BEGIN { OFS="\t"; beta_col=0; se_col=0 }

NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == beta_name) beta_col=i
    if ($i == se_name)   se_col=i
  }
  if (beta_col==0 || se_col==0) {
    print "ERROR: couldn't find beta/se columns. beta_name=" beta_name ", se_name=" se_name #> "/dev/stderr"
    exit 2
  }
  print $0, "Z"
  next
}

{
  beta = $beta_col
  se   = $se_col
  z = (beta=="" || se=="" || beta=="NA" || se=="NA" || se==0) ? "NA" : beta/se
  print $0, z
}
