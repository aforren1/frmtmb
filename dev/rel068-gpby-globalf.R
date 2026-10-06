# test-gp-by.R with a global `f` defined, as the release runner had one
# before it ran inside local() (the release review, m6): the refusal
# test must not depend on the global environment.
#   Rscript dev/rel068-gpby-globalf.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
f <- "C:/a/path/that/answers/for/f"
g <- "another global"
suppressMessages({library(testthat); library(frmtmb)})
r <- as.data.frame(test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-gp-by.R",
  package = "frmtmb", env = test_env("frmtmb"), reporter = "silent"))
cat("global f:", exists("f", globalenv()), " pass", sum(r$passed), "fail",
    sum(r$failed), "err", sum(r$error), "skip", sum(r$skipped), "\n")
