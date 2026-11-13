# Lane gpstop: one test file in one R process. Args: arm (lane | base),
# package, test file path. Prints the lib: line and a RESULT line.
a <- commandArgs(TRUE)
arm <- a[1]; p <- a[2]; f <- a[3]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpstop-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/stan-cache")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = p)), p,
    as.character(packageVersion(p)), "| frmtmb",
    dirname(system.file(package = "frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
