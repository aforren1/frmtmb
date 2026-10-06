# One line per firing of the standard-error warning in a
# dev/nanse-diag.sh directory: the test file, the lost parameters with
# their reason, and what base's sdreport() gave for the same fit.
#   Rscript dev/nanse-fire-sum.R <dir>
d <- commandArgs(TRUE)[1]
fs <- list.files(d, pattern = "fire$", full.names = TRUE)
n_all <- 0L
tab <- list()
for (f in fs) {
  x <- readLines(f)
  i <- grep("^Standard errors are not available", x)
  for (k in i) {
    msg <- sub(" The other standard errors.*", "", x[k])
    what <- sub("^Standard errors are not available for ", "", msg)
    reason <- if (grepl("a bound holds", what)) "bound" else
      if (grepl("curves downward", what)) "concave" else
      if (grepl("not finite in", what)) "nonfinite" else "flat"
    pars <- sub(":.*$", "", sub("^[0-9]+ of [0-9]+ parameters[.] ", "", what))
    base <- sub(" [|] random.*", "", sub("^  base: ", "", x[k + 2L]))
    tab[[length(tab) + 1L]] <- data.frame(
      file = sub("[.]fire$", "", basename(f)),
      lost = paste0(sub(" parameters.*", "", what), " (", reason, ")"),
      base = base, stringsAsFactors = FALSE)
    n_all <- n_all + 1L
  }
}
X <- do.call(rbind, tab)
cat("firings:", n_all, " base sdreport() non-finite on:",
    sum(!grepl("^0 of", X$base)), "\n")
options(width = 200)
why <- sub("[)]$", "", sub(".*[(]", "", X$lost))
cat("by reason:", paste(names(table(why)), table(why), collapse = "; "),
    "\n")
print(X, row.names = FALSE, right = FALSE)
