# Review of lane famend: one test file of core on the lane build, under
# whichever R runs this script (the OpenBLAS copies of
# dev/ciharden-openblas.sh).
#   Rscript dev/reviews/famend-review-runtest.R <file>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("R", R.home(), "\n")
cat("BLAS", extSoftVersion()[["BLAS"]], "LAPACK", La_library(), "\n")
cat("lib", dirname(find.package("frmtmb")), "\n")
f <- file.path("C:/Users/adf44/source/r/frmtmb-wt-famend/tests/testthat", a[1])
res <- testthat::test_file(f, package = "frmtmb",
                           env = testthat::test_env("frmtmb"),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("FILE %s pass %d fail %d err %d skip %d warn %d\n", a[1],
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped),
            sum(df$warning)))
for (i in which(df$failed > 0 | df$error | df$warning > 0)) {
  cat("  not ok:", df$test[i], " failed", df$failed[i], "error", df$error[i],
      "warn", df$warning[i], "\n")
}
for (r in res) for (e in r$results) {
  if (!inherits(e, "expectation_success") && !inherits(e, "expectation_skip")) {
    cat("  --", r$test, ":", substr(conditionMessage(e), 1, 400), "\n")
  }
}
