# One test file, one R process. Counts come from the result data frame's
# own columns, including `error` and `skipped`, because a runner that sums
# only `failed` prints a clean line for a file that aborted halfway.
#
#   Rscript dev/gradcheck-runtests.R lane tests/testthat/test-x.R ...

args <- commandArgs(trailingOnly = TRUE)
which_lib <- args[[1L]]
files <- args[-1L]
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
Sys.setenv(NOT_CRAN = "true")
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n")
for (f in files) {
  res <- tryCatch(
    as.data.frame(testthat::test_file(
      f, package = "frmtmb", env = testthat::test_env("frmtmb"),
      reporter = "silent")),
    error = function(e) {
      cat(sprintf("%-42s ABORTED %s\n", basename(f),
                  conditionMessage(e)))
      NULL
    })
  if (is.null(res)) next
  cat(sprintf("%-42s pass=%d fail=%d error=%d skip=%d\n", basename(f),
              sum(res$passed), sum(res$failed), sum(res$error),
              sum(res$skipped)))
  bad <- res[res$failed > 0 | res$error > 0, , drop = FALSE]
  if (nrow(bad)) {
    for (i in seq_len(nrow(bad))) {
      cat("   x ", bad$test[i], "\n", sep = "")
    }
  }
  utils::flush.console()
}
