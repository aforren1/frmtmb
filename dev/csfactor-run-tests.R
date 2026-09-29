# One test file per R process, counts summed from the result frame.
#   Rscript dev/csfactor-run-tests.R <lib> <file> [<file> ...]
# The library comes first so a "before" run against rellib-r3 and an
# "after" run against the lane library differ only in that argument.
args <- commandArgs(TRUE)
LIB <- args[[1L]]
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
WT <- "C:/Users/adf44/source/r/frmtmb-wt-csfactor"
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
for (nm in args[-1L]) {
  f <- file.path(WT, "tests", "testthat", nm)
  res <- as.data.frame(testthat::test_file(
    f, package = "frmtmb", env = testthat::test_env("frmtmb"),
    reporter = "silent"))
  cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n", nm,
              sum(res$passed), sum(res$failed), sum(res$error),
              sum(res$skipped), sum(res$warning)))
  bad <- res[res$failed > 0 | res$error, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat("  BAD:", bad$test[i], "\n")
    for (r in bad$result[[i]]) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        cat("    ", gsub("\n", "\n     ", conditionMessage(r)), "\n")
      }
    }
  }
  sk <- res[res$skipped, , drop = FALSE]
  for (i in seq_len(nrow(sk))) cat("  SKIP:", sk$test[i], "\n")
}
cat("done\n")
