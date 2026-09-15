source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# The review measured the SWAP build at 22 of 26 in mode U; this script
# measured 23 of 26. Rather than argue which is right, print the names,
# so the one-name difference can be attributed.
sub1("dev/generics-propagate.R",
paste0("cat(sprintf(\"re-exported from core      lost %d of %d\\n\",\n",
       "            length(a), length(reexported)))\n"),
paste0("cat(sprintf(\"re-exported from core      lost %d of %d\\n\",\n",
       "            length(a), length(reexported)))\n",
       "# which of the lost ones are answered by SOMEBODY'S default rather\n",
       "# than by nothing: a count that only looks for \"no method at all\"\n",
       "# will not count these, which is the method difference that makes\n",
       "# the review's 22 and this script's 23 disagree by one\n",
       "dflt <- a[vapply(a, reach, NA, \"default\")]\n",
       "cat(sprintf(\"  of which answered by a .default  %d  %s\\n\",\n",
       "            length(dflt), paste(dflt, collapse = \",\")))\n"))
cat("DONE\n")
