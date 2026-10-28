# Reviewer: one test file in one process.
#   Rscript dev/review-priordraws-runtest.R <pkg> <test file> <lane|base>
a <- commandArgs(TRUE)
pkg <- a[1]
f <- a[2]
arm <- if (length(a) > 2) a[3] else "lane"
libs <- c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-priordraws-lib",
          if (identical(arm, "rv")) c("C:/Users/adf44/source/r/wt-priordraws-rv-lib",
                                      "C:/Users/adf44/source/r/wt-priordraws-lib"),
          "C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), "|", find.package(pkg), "\n")
ws <- character(0)
r <- withCallingHandlers(
  testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                      reporter = "summary"),
  warning = function(w) ws <<- c(ws, conditionMessage(w)))
df <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d escaped=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning), length(ws)))
if (length(ws)) cat("ESCAPED:", unique(ws), sep = "\n  ")
