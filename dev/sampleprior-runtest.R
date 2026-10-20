# Lane sampleprior: run one test file of one package against a library,
# one file per process, as dev/lane-rules.md asks.
#   Rscript dev/sampleprior-runtest.R <lane|r9> <pkg> <file>
a <- commandArgs(TRUE)
arm <- a[1]
pkg <- a[2]
file <- a[3]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/stan-cache")
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package(pkg), " frmtmb:", find.package("frmtmb"), "\n")
r <- testthat::test_file(file, package = pkg,
                         env = testthat::test_env(pkg),
                         reporter = testthat::ProgressReporter$new(
                           show_praise = FALSE))
d <- as.data.frame(r)
cat(sprintf("RESULT %s: tests %d pass %d fail %d error %d skip %d warn %d\n",
            basename(file), nrow(d), sum(d$passed), sum(d$failed),
            sum(d$error), sum(d$skipped), sum(d$warning)))
