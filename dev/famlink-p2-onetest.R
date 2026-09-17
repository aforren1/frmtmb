# Run one test file against a library. Usage:
#   Rscript dev/famlink-p2-onetest.R <lib> <file> [filter]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat)
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
r <- as.data.frame(test_file(a[2], reporter = "progress",
                             package = "frmtmb", load_package = "installed"))
cat(sprintf("FAMLINK %s tests %d fail %d error %d skip %d\n", basename(a[2]),
            nrow(r), sum(r$failed), sum(r$error), sum(r$skipped)))
