# Lane nlstop2: run one test file of a package on the arm's build, as
# lane-env-r083.md says (library() first, test_env, silent reporter), and
# print each failure or error and a RESULT line with pass, fail, error,
# skip and warning counts.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-testfile.R <file>
#   [package]
a <- commandArgs(TRUE)
f <- a[1]
p <- if (length(a) > 1) a[2] else "frmtmb"
arm <- Sys.getenv("NLSTOP2_ARM", "after")
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-nlstop2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-nlstop2/dev/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "arm:", arm, "\n")
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "silent")
for (t in res) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning"))) {
      cat("FAILED [", class(e)[1], "]", t$test, ":",
          substr(gsub("\\s+", " ", conditionMessage(e)), 1, 400), "\n")
    }
  }
}
r <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d secs=%.1f\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning),
            proc.time()[["elapsed"]] - t0))
