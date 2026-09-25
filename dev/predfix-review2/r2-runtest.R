.libPaths(c("C:/Users/adf44/source/r/predfix-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
args <- commandArgs(TRUE)
f <- args[1]; pkg <- args[2]
suppressPackageStartupMessages({library(testthat); library(pkg, character.only = TRUE)})
cat("frmtmb from", find.package("frmtmb"), "\n")
dir <- dirname(f)
before <- file.exists(file.path(dir, "Rplots.pdf"))
r <- as.data.frame(testthat::test_file(f, reporter = "silent",
                                      package = pkg, load_package = "installed"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", basename(f),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
cat("Rplots.pdf before:", before, "after:", file.exists(file.path(dir, "Rplots.pdf")), "\n")
