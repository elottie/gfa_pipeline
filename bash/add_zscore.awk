# add_zscore.awk
BEGIN { OFS="\t"; beta_col=-1; se_col=-1 }
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="beta") beta_col=i
        if($i=="se") se_col=i
    }
    print $0 "\tz"
    next
}
NR>1 {
    beta=$beta_col
    se=$se_col
    z = (beta=="" || se=="" || beta=="NA" || se=="NA" || se==0) ? "NA" : beta/se
    print $0, z
}
