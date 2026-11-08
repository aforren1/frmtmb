# Review of lane mfxse: run one test file of core on one library set, the
# way the release harness does (package attached, namespace as test env).
#   Rscript dev/review-mfxse-run-test.R <file> [lib]
#   lib: "lane" (default, wt-mfxse-lib then rellib-r18), "before"
#        (rellib-r18 alone), or a path put before rellib-r18 (a mutant).
a <- commandArgs(TRUE)
f <- a[1]
arm <- if (length(a) > 1) a[2] else "lane"
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(switch(arm, lane = c(lane, base, user), before = c(base, user),
                 c(arm, base, user)))
Sys.setenv(NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "marginaleffects", format(packageVersion("marginaleffects")), "\n")
res <- as.data.frame(testthat::test_file(f, package = "frmtmb",
                                         env = testthat::test_env("frmtmb"),
                                         reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(res$passed), sum(res$failed),
            sum(res$error), sum(res$skipped), sum(res$warning)))
