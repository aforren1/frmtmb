# The smallfix review's bitwise sample, base (rellib-r12) against merged
# (rellib-r13), each twice (dev/rel075-eq-smallfix.bat). Adapted from
# the review's r1-compare.R: identical() on each saved output, and for a
# difference its largest relative difference, or the first differing
# line of printed text. A difference between base1 and base2 is the
# process's own (rtmb-pitfalls item 21, or an environment in a saved
# object); a difference between base and merged that the base arms do
# not show is a change of the build, to be accounted for by a lane.
#   Rscript dev/rel075-eq-smallfix.R
O <- "C:/Users/adf44/source/r/frmtmb-wt-rel075/dev/rel075-log/eq/smallfix/"
r <- lapply(c(base1 = "base1", base2 = "base2", merged1 = "merged1",
              merged2 = "merged2"), function(s) {
  readRDS(paste0(O, s, ".rds"))
})
cmp <- function(x, y) {
  if (identical(x, y)) return("identical")
  ux <- unlist(x); uy <- unlist(y)
  if (is.numeric(ux) && is.numeric(uy) && length(ux) == length(uy)) {
    m <- max(abs(ux - uy) / pmax(1, abs(ux)), na.rm = TRUE)
    return(sprintf("DIFFER (max rel %.3g)", m))
  }
  if (is.character(ux) && is.character(uy)) {
    k <- which(ux[seq_len(min(length(ux), length(uy)))] !=
                 uy[seq_len(min(length(ux), length(uy)))])[1]
    if (is.na(k)) return(sprintf("DIFFER (lengths %d, %d)", length(ux),
                                 length(uy)))
    return(sprintf("DIFFER [%s] -> [%s]", substr(ux[k], 1, 90),
                   substr(uy[k], 1, 90)))
  }
  "DIFFER"
}
pairs <- list(c("base1", "base2"), c("merged1", "merged2"),
              c("base1", "merged1"), c("base2", "merged2"))
nout <- 0L; nsame <- 0L
for (nm in names(r$base1)) {
  for (o in names(r$base1[[nm]])) {
    nout <- nout + 1L
    res <- vapply(pairs, function(p) cmp(r[[p[1]]][[nm]][[o]],
                                         r[[p[2]]][[nm]][[o]]), "")
    if (all(res == "identical")) { nsame <- nsame + 1L; next }
    cat(sprintf("%-11s %-9s %s\n", nm, o,
                paste(sapply(pairs, paste, collapse = "~"), res,
                      sep = ": ", collapse = " | ")))
  }
}
cat("outputs", nout, "| identical in all four pairs", nsame, "\n")
