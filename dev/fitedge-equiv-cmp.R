# The remeasured pairs (dev/fitedge-equiv-out) against the record
# (dev/equiv-data/results.tsv): verdict, level, and every field that moved.
o <- read.delim("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/equiv-data/results.tsv")
for (f in list.files("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-equiv-out",
                     "[.]tsv$", full.names = TRUE)) {
  n <- read.delim(f)
  r <- o[o$id == n$id, ]
  k <- intersect(names(r), names(n))
  moved <- k[vapply(k, function(j) {
    !identical(format(r[[j]]), format(n[[j]]))
  }, NA)]
  cat(sprintf("EQ %s verdict %s -> %s, level %s -> %s\n", n$id, r$verdict,
              n$verdict, r$level, n$level))
  for (j in moved) cat("EQ   ", j, ":", format(r[[j]]), "->", format(n[[j]]), "\n")
}
