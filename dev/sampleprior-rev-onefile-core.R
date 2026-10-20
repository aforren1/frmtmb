# Reviewer of lane sampleprior: one test file of core frmtmb per
# process, gated, against the lane library.
#   Rscript dev/sampleprior-rev-onefile.R <file>
f <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), " frmtmb:",
    find.package("frmtmb"), "\n")
d <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/tests/testthat"
r <- testthat::test_file(file.path(d, f), package = "frmtmb",
                         env = testthat::test_env("frmtmb"),
                         reporter = "summary")
x <- as.data.frame(r)
cat(sprintf("RESULT %s: tests %d pass %d fail %d error %d skip %d warn %d\n",
            f, nrow(x), sum(x$passed), sum(x$failed), sum(x$error),
            sum(x$skipped), sum(x$warning)))
