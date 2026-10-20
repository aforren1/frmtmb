# One test file, one R process, on the lane library or the base.
#   Rscript dev/catmulti-test1.R <arm: lane|base> <package> <test file>
# Prints a RESULT line with pass/fail/error/skip/warn and the library.
a <- commandArgs(TRUE)
arm <- a[1]; p <- a[2]; f <- a[3]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(lib, base, user) else c(base, user))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(NOT_CRAN = "true")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache")
}
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "| pkg",
    dirname(find.package(p)), "\n")
res <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                           reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d | core %s\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning),
            dirname(find.package("frmtmb"))))
bad <- df[df$failed > 0 | df$error, c("test", "failed", "error")]
if (nrow(bad)) {
  cat("NOT PASSING (every one, uncapped):\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %s | failed=%d error=%s\n", bad$test[i], bad$failed[i],
                bad$error[i]))
  }
}
