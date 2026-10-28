# One test file in one process, against a chosen library stack.
#   Rscript dev/priordraws-runtest.R <pkg> <test file> <lane|base>
# "lane" puts the lane's library first; "base" reads rellib-r14 alone
# (the before arm). Prints one RESULT line with the counts.
a <- commandArgs(TRUE)
pkg <- a[1]
f <- a[2]
arm <- if (length(a) > 2) a[3] else "lane"
libs <- c(if (identical(arm, "lane")) "C:/Users/adf44/source/r/wt-priordraws-lib",
          "C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/stan-cache")
}
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), "|", find.package(pkg), "\n")
r <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                         reporter = "summary")
df <- as.data.frame(r)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
