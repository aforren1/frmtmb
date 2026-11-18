# Lane squant: run one test file against the lane library.
# Run: Rscript dev/squant-test1.R <package> <test file path>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- a[1]
suppressMessages({library(testthat); library(p, character.only = TRUE)})
cat("lib:", find.package(p), format(packageVersion(p)), "frmtmb at",
    find.package("frmtmb"), "\n")
t0 <- proc.time()[["elapsed"]]
r <- as.data.frame(testthat::test_file(a[2], package = p,
                                       env = testthat::test_env(p),
                                       reporter = "silent"))
fails <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(fails)) print(fails)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d %.1fs\n",
            basename(a[2]), sum(r$passed), sum(r$failed), sum(r$error),
            sum(r$skipped), sum(r$warning), proc.time()[["elapsed"]] - t0))
