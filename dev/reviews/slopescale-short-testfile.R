# Short check of lane slopescale's fix round: one test file on the lane
# library. Usage: Rscript slopescale-short-testfile.R <test file name>
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- "frmtmb"
library(p, character.only = TRUE)
cat("lib:", find.package(p), "\n")
f <- file.path("C:/Users/adf44/source/r/frmtmb-wt-slopescale/tests/testthat",
               args[1])
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
cat("RESULT", args[1], "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning),
    "\n")
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("BAD:", bad, sep = "\n  ")
