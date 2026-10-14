# Lane illse fix round 1: one test file against the lane library or,
# with ILLSE_ARM=base, rellib-r8 alone; the counts per test.
#   Rscript dev/illse-fx-tbytest.R <test file path>
arm <- Sys.getenv("ILLSE_ARM", "lane")
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-illse-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
res <- as.data.frame(test_file(commandArgs(TRUE)[1], reporter = "silent",
                               package = "frmtmb", load_package = "none"))
res$test <- substr(res$test, 1, 50)
print(res[, c("test", "nb", "passed", "failed", "error", "warning")],
      right = FALSE, row.names = FALSE, width = 200)
