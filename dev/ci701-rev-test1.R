# Review of lane ci701: run one core test file on the lane build.
#   Rscript dev/ci701-rev-test1.R <file> [lane|base]
a <- commandArgs(TRUE)
f <- a[1]
arm <- if (length(a) > 1) a[2] else "lane"
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "base") c(base, ul) else
  c("C:/Users/adf44/source/r/wt-ci701-lib", base, ul))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
r <- as.data.frame(testthat::test_file(
  file.path("C:/Users/adf44/source/r/frmtmb-wt-ci701/tests/testthat", f),
  package = "frmtmb", env = testthat::test_env("frmtmb"),
  reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n", f,
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
            sum(r$warning)))
bad <- r[r$failed > 0 | r$error | r$warning > 0, c("test", "failed",
                                                     "warning")]
if (nrow(bad)) print(bad)
