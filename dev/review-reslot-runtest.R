# Reviewer's one-file test runner. Usage:
# Rscript dev/review-reslot-runtest.R before|after <path> <package>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
f <- args[2]
p <- args[3]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "before")) base else
  c("C:/Users/adf44/source/r/wt-reslot-lib", base))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/stan-cache")
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "|", p,
    dirname(find.package(p)), "\n")
library(testthat)
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "summary"))
cat(sprintf("RESULT %s %s pass=%d fail=%d error=%d skip=%d warn=%d\n", arm,
            basename(f), sum(r$nb) - sum(r$failed) - sum(r$skipped),
            sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
