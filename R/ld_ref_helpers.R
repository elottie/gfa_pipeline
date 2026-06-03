awk_dedup_ld_ref <- '
  {count[$1]++
  if (count[$1] == 1) {
     order[++n] = $1
   line[$1] = $0
      }
    }
  END {
	      for (i = 1; i <= n; i++) {
		            snp = order[i]
        if (count[snp] == 1) {
		        print line[snp]
	      }
	    }
    }
'

build_concat_ld_ref_cmd <- function(ld_files,
				                                        snps_in_ref_file,
									                                    awk_get_snp_and_l2,
									                                    awk_dedup = awk_dedup_ld_ref) {
	  sprintf(
		      "{ printf 'snp\\tl2\\n'; for f in %s; do zcat -- \"$f\" | awk -F%s -v file=\"$f\" %s; done | awk -F%s %s; } > %s",
		          paste(shQuote(ld_files), collapse = " "),
		          shQuote("\t"),
			      shQuote(awk_get_snp_and_l2),
			      shQuote("\t"),
			          shQuote(awk_dedup),
			          shQuote(snps_in_ref_file)
				    )
}
