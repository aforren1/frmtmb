# The failures test-grad-verdict.R produces on the 0.64.0 reference
# build, expectation by expectation. A test that pins a defect is
# worthless unless it has been SEEN TO FAIL, and an error because a
# symbol does not exist is the weak form, so the message of each one is
# recorded here.
#
#   Rscript dev/gradcheck-10-seenfail.R base|lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
Sys.setenv(NOT_CRAN = "true")
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

res <- as.data.frame(testthat::test_file(
  "tests/testthat/test-grad-verdict.R", package = "frmtmb",
  env = testthat::test_env("frmtmb"), reporter = "silent"))
cat(sprintf("pass=%d fail=%d error=%d skip=%d\n", sum(res$passed),
            sum(res$failed), sum(res$error), sum(res$skipped)))
for (i in seq_len(nrow(res))) {
  if (res$failed[i] == 0 && res$error[i] == 0) next
  cat("\n== ", res$test[i], "\n", sep = "")
  for (r in res$result[[i]]) {
    cls <- class(r)[1L]
    if (!cls %in% c("expectation_failure", "expectation_error")) next
    cat("   [", cls, "] ",
        gsub("[\r\n]+", " | ", substr(conditionMessage(r), 1, 240)),
        "\n", sep = "")
  }
}
