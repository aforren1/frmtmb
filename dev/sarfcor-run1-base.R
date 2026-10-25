# Run one test file against the BASE library alone (rellib-r12), to see
# a test fail before the change: Rscript dev/sarfcor-run1-base.R <pkg> <file>
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- args[1]
suppressMessages(library(p, character.only = TRUE))
cat("lib:", p, find.package(p), format(packageVersion(p)), "| frmtmb",
    find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(args[2], package = p,
                                       env = testthat::test_env(p),
                                       reporter = "summary"))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n",
            basename(args[2]), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning)))
