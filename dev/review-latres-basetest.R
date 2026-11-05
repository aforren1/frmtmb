# Attack 10: the lane's test files against the before-build, rellib-r16. Which
# failures are the behavior, and which are only a missing symbol?
Sys.setenv(REV_ARM = "base")
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages(library(testthat))
print(packageVersion("frmtmb"))
res <- test_file("tests/testthat/test-latres.R", package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("tests", sum(df$nb), "passed", sum(df$passed), "failed", sum(df$failed), "errors", sum(df$error), "\n")
for (r in res) for (e in r$results) {
  if (inherits(e, "expectation_success")) next
  cat(sprintf("[%s] %s\n", substr(r$test, 1, 55), substr(gsub("\n", " ", conditionMessage(e)), 1, 150)))
}
