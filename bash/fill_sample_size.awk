# fill_sample_size.awk
# expects -v ss_name="$sample_size" -v pub_ss_val="$pub_sample_size"

BEGIN {
    FS = OFS = "\t"
    ss_col = -1
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == ss_name) ss_col = i
    }

    if (ss_col == -1) {
        print "ERROR: Could not find sample size column: " ss_name #> "/dev/stderr"
        exit 1
    }

    if (pub_ss_val == "" || pub_ss_val == "NA"){
	print "ERROR: Could not find published sample size value: " pub_ss_val
	exit 1
    #} else {
    #    print "Found published sample size value, filling in sample_size NAs with it: " pub_ss_val #> "/dev/stderr"
    }

    print
    next
}

{
    if ($ss_col == "" || $ss_col == "NA") {
        $ss_col = pub_ss_val
    }

    print
}
