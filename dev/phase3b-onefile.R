# Run one test file against a named library first on the path.
# Usage: Rscript dev/phase3b-onefile.R <lib> <pkg> <test file>
a <- commandArgs(TRUE)
.libPaths(unique(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
                   "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/source/r/pinlib",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(testthat); library(a[2], character.only = TRUE)})
cat("from", find.package(a[2]), "\n")
r <- as.data.frame(test_file(a[3], reporter = "progress", package = a[2],
                             load_package = "installed"))
cat("RESULT", basename(a[3]), "blocks", nrow(r), "expectations", sum(r$nb),
    "pass", sum(r$passed), "fail", sum(r$failed), "err", sum(r$error), "\n")
