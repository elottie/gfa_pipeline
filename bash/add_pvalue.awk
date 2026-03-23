# add_pvalue.awk
BEGIN { OFS="\t"; z_col=0 }

NR==1 {
  for (i=1;i<=NF;i++) if ($i=="Z" || $i=="z") z_col=i
  print $0, "p"
  next
}

{
  z = (z_col? $z_col : "NA")
  if (z=="" || z=="NA") p="NA";
  else {
    az = (z<0 ? -z : z)
    # two-sided p = erfc(|z|/sqrt(2))
    p = erfc(az / sqrt(2))
  }
  print $0, p
}
