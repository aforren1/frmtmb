# Reviewer of lane optrap: one test file with lane optrap's response units
# taken out in this process (no unit: predictor_units() empty), to tell
# an interplay with the units from anything else on a merged build.
#   Rscript dev/optrap-rev-ablate-tf.R <lib> <ablate: units|none> <file>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
if (a[2] == "units") {
  utils::assignInNamespace("predictor_units", function(frame) numeric(0),
                           "frmtmb")
}
cat("ablate", a[2], "lib", dirname(find.package("frmtmb")), "\n")
res <- test_file(a[3], package = "frmtmb",
                 env = testthat::test_env("frmtmb"), reporter = "silent")
r <- as.data.frame(res)
for (t in res) for (x in t$results) {
  if (inherits(x, c("expectation_failure", "expectation_error"))) {
    cat("  FAIL [", t$test, "]", substr(conditionMessage(x), 1, 200), "\n")
  }
}
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d warn=%d\n",
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
