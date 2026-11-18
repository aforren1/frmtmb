# Review of lane quadesc: run tests/testthat/test-quad-escape-record.R of
# the worktree against <lib|base>, silent reporter, results listed.
#   Rscript dev/reviews/quadesc-review-testfile.R <lib|base>
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- "frmtmb"
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package(p)), "\n")
t0 <- proc.time()[["elapsed"]]
r <- testthat::test_file("tests/testthat/test-quad-escape-record.R",
                         package = p, env = testthat::test_env(p),
                         reporter = "silent")
df <- as.data.frame(r)
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")])
for (res in r) {
  for (x in res$results) {
    if (!inherits(x, "expectation_success")) {
      cat("NOT OK [", res$test, "]", conditionMessage(x), "\n")
    }
  }
}
cat(sprintf("elapsed %.1f s\n", proc.time()[["elapsed"]] - t0))
