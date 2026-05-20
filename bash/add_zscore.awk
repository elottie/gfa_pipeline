# add_zscore.awk
# expects -v beta_name="..." -v se_name="..."

# real awk ppl like to initialize to 0 bc awk starts with 1.  but I think negative 1 is even clearer & def gives error
BEGIN {
    FS=OFS="\t"
    beta_col=-1
    se_col=-1
}

NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == beta_name) beta_col=i
    if ($i == se_name)   se_col=i
  }
  if (beta_col==-1 || se_col==-1) {
    print "ERROR: couldn't find beta/se columns. beta_name=" beta_name ", se_name=" se_name #> "/dev/stderr"
    print "DEBUG NF=" NF #> "/dev/stderr"
        for (i = 1; i <= NF; i++) {
            print "DEBUG field " i " = [" $i "]" #> "/dev/stderr"
        }
    exit 1
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
