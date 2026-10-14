# Reviewer, lane dblbar: run one test file against the lane or the base.
# Usage: Rscript dev/dblbar-rev-runtest.R base|lane <pkg> <file>
a <- commandArgs(TRUE)
arm <- a[1]; p <- a[2]; f <- a[3]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-dblbar/dev/stan-cache")
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), find.package(p), "\n")
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d skip=%d warn=%d\n", arm,
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning)))
bad <- r[r$failed > 0 | r$error, ]
if (nrow(bad)) print(bad[, c("test", "failed", "error")])
