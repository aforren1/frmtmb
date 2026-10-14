# Lane illse: run one core test file against the worktree as it is
# (pkgload::load_all(), pure R), for iterating while the lane library is
# in use by a scan. Not a substitute for the installed runs.
#   Rscript dev/illse-devtest.R <test file>
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-illse",
                                   quiet = TRUE))
suppressMessages(library(testthat))
res <- as.data.frame(test_file(commandArgs(TRUE)[1], reporter = "summary",
                               package = "frmtmb", load_package = "none"))
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d warn=%d\n",
            sum(res$passed), sum(res$failed), sum(res$error),
            sum(res$skipped), sum(res$warning)))
