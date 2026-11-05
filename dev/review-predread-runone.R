# One test file, one process, with a given core ahead of the lane's library.
#   Rscript dev/review-predread-runone.R <package> <file> <core lib>
a <- commandArgs(TRUE)
.libPaths(c(a[3], "C:/Users/adf44/source/r/wt-predread-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
p <- a[1]
suppressMessages(library(p, character.only = TRUE))
cat("lib:", p, dirname(find.package(p)), "; frmtmb", dirname(find.package("frmtmb")), "\n")
r <- as.data.frame(testthat::test_file(a[2], package = p, env = testthat::test_env(p),
                                       reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n", basename(a[2]),
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning)))
bad <- r[r$failed > 0 | r$error, c("test")]
if (length(bad)) cat("FAILED:", paste(unique(bad), collapse = " || "), "\n")
