# test-conditions.R's "immediate." block runs Rscript in a subprocess,
# and Rscript has been seen to segfault on exit on this machine, which
# fails that one block. Run the file on a chosen library, attached as
# tests/testthat.R attaches it, and print one line per run.
#
#   Rscript dev/rel066-conditions-flake.R <library> <label>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
f <- "tests/testthat/test-conditions.R"
r <- as.data.frame(test_file(f, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent"))
bad <- r$test[r$failed > 0 | r$error]
cat("FLAKE ", a[2], " frmtmb ", as.character(packageVersion("frmtmb")),
    " pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
    sum(r$error), if (length(bad)) paste0(" failing: ", bad), "\n",
    sep = "")
