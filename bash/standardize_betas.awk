#awk -v beta_name="$beta_hat" -v effect_or_flag="$effect_or" -f standardize_betas.awk

# real awk ppl like to initialize to 0 bc awk starts with 1.  but I think negative 1 is even clearer & def gives error

BEGIN {
    FS = OFS = "\t"
    beta_col = -1
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == beta_name) beta_col = i
    }

    if (beta_col == -1) {
        print "ERROR: Could not find beta hat column: " beta_name
        exit 1
    }

    effect_or_flag_upper=toupper(effect_or_flag)
    if ( effect_or_flag_upper != "TRUE" && effect_or_flag_upper != "T" && effect_or_flag_upper != "FALSE" && effect_or_flag_upper != "F" ) {
        print "ERROR: Invalid effect OR flag: " effect_or_flag
        exit 1
    #} else {
    #    print "Found effect OR flag, ensuring betas are log(OR) with it: " effect_or_flag #> "/dev/stderr"
    }

    print
    next
}

{
    if (effect_or_flag_upper == "TRUE" || effect_or_flag_upper == "T") {
        # beta is OR, so convert to log(OR)
	if (($beta_col + 0) > 0) {
            $beta_col = log($beta_col + 0)
        } else {
            $beta_col = "NA"
        }
    }
    print
}
