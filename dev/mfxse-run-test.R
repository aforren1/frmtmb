# Run one test file of one package on one arm, the way the release harness
# does: the package attached, its namespace as the test environment.
#   MFX_ARM=before|after Rscript dev/mfxse-run-test.R <pkg> <file>
a <- commandArgs(TRUE)
pkg <- a[1]
f <- a[2]
arm <- Sys.getenv("MFX_ARM", "after")
lib <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "before") c(base, user) else c(lib, base, user))
Sys.setenv(NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/stan-cache",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(testthat)
  library(pkg, character.only = TRUE)
})
cat("lib:", find.package("frmtmb"), " pkg:", find.package(pkg), "\n")
res <- as.data.frame(testthat::test_file(f, package = pkg,
                                         env = testthat::test_env(pkg),
                                         reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(res$passed), sum(res$failed),
            sum(res$error), sum(res$skipped), sum(res$warning)))
