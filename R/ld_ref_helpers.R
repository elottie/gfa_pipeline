# make awks, sorts, and joins consistent across users
Sys.setenv(LC_ALL = "C")

# NR > 1 so we do not copy header
# assuming tab-sep input file
awk_get_snp_and_l2 <- '
  BEGIN { OFS = FS }
  NR == 1 {
    snp = -1
    l2 = -1

    for (i = 1; i <= NF; i++) {
      if ($i == "SNP") snp = i
      if ($i == "L2")  l2 = i
    }

    if (snp == -1 || l2 == -1) {
      print "Missing SNP or L2 column in " file #> "/dev/stderr"
      exit 1
    }

    next
  }

  {
    print $snp, $l2
  }
'

# this is seeming more complicated.  it is really just removing all instances of a snps that is seen more than once (vs doing unique)
# don't need to set ofs bc doing print line
# don't need to sort anymore bc not doing the join/sort approach
awk_dedup_ld_ref <- '
  {
    count[$1]++
  
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

build_concat_ld_ref <- function(ld_files,
		                snps_in_ref_file,
			        delim="\t") {
  sprintf(
    "{ printf 'snp\\tl2\\n'; for f in %s; do zcat -- \"$f\" | awk -F%s -v file=\"$f\" %s; done | awk -F%s %s; } > %s",
    paste(shQuote(ld_files), collapse = " "),
    shQuote(delim),
    shQuote(awk_get_snp_and_l2),
    shQuote(delim),
    shQuote(awk_dedup_ld_ref),
    shQuote(snps_in_ref_file)
  )
}
