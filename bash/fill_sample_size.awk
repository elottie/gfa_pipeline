# fill_sample_size.awk
BEGIN { OFS="\t"; ss_col=-1 }
NR==1 {
    for(i=1;i<=NF;i++)
        if($i=="sample_size") ss_col=i
    print
}
NR>1 {
    if($ss_col == "" || $ss_col == "NA") $ss_col=pub_ss
    print
}
