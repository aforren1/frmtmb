# One test file per process on the lane library. Args: package, file.
args <- commandArgs(TRUE)
p <- args[[1L]]
f <- args[[2L]]
.libPaths(c("C:/Users/adf44/source/r/wt-refits-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
t0 <- proc.time()[["elapsed"]]
r <- as.data.frame(testthat::test_file(f, package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
el <- proc.time()[["elapsed"]] - t0
bad <- r[r$failed > 0 | r$error, "test"]
cat(sprintf("%-32s pass=%d fail=%d error=%d skip=%d warn=%d tests=%d %.0fs lib:%s %s\n",
            basename(f), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning), nrow(r), el,
            dirname(find.package("frmtmb")),
            dirname(find.package(p))))
if (length(bad)) cat("  FAILED:", paste(bad, collapse = " | "), "\n")
