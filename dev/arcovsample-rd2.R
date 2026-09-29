# Lane wt-arcovsample, punch round 1: RENDER the two corrected topics
# and read the corrected sentence out of the RENDERED text, not out of
# the source. The blocking defect was in a rendered sentence, so the
# check has to be on the rendering (dev/lane-rules.md).
#
#   Rscript dev/arcovsample-rd2.R > dev/arcovsample-log/rd2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-arcovsample-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

show <- function(f, key, before = 2L, after = 8L) {
  o <- tempfile(fileext = ".txt")
  tools::Rd2txt(f, out = o)
  t <- readLines(o, warn = FALSE)
  cat("\n==== ", basename(f), " (", length(t), " rendered lines)\n",
      sep = "")
  i <- grep(key, t, fixed = TRUE)
  if (!length(i)) {
    cat("  NOT FOUND in the rendering: ", key, "\n", sep = "")
    return(invisible(NULL))
  }
  k <- sort(unique(unlist(lapply(i, function(x) {
    seq(max(1L, x - before), min(length(t), x + after))
  }))))
  cat(paste(t[k], collapse = "\n"), "\n")
  # the refuted wording must be GONE from the rendering
  bad <- grep("first `max(p, q)` rows", t, fixed = TRUE)
  bad2 <- grep("first max(p, q) rows", t, fixed = TRUE)
  cat("  refuted wording still present: ",
      length(bad) + length(bad2) > 0L, "\n", sep = "")
  cat("  '%' surviving into the rendering: ",
      grepl("%", paste(t, collapse = ""), fixed = TRUE), "\n", sep = "")
  invisible(NULL)
}

show("extensions/frmtmb.sample/man/sample-log_lik.Rd", "lag")
show("man/frmtmb-autocor.Rd", "CONDITIONAL likelihood")
show("man/frmtmb-autocor.Rd", "POINTWISE")
