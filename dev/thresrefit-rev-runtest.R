## Reviewer's own runner: ONE test file per R process. Counts come from
## the per-expectation data frame, not from the reporter's summary, which
## caps at ten failures. It also prints the number of test blocks and any
## file that reported nothing, so an aborted file cannot pass as clean.
## Usage: FRMTMB_LIB=[lane|base] Rscript dev/thresrefit-rev-runtest.R <f>
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/wt-thresrefit-lib")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
f <- commandArgs(trailingOnly = TRUE)[1L]
root <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit"
library(frmtmb)
library(testthat)
setwd(file.path(root, "tests/testthat"))
res <- test_file(file.path(root, "tests/testthat", f), package = "frmtmb",
                 env = test_env("frmtmb"), reporter = "silent")
df <- as.data.frame(res)
cat("REVRESULT ", f, " arm=", arm,
    " blocks=", nrow(df),
    " pass=", sum(df$passed),
    " fail=", sum(df$failed),
    " err=", sum(df$error),
    " warn=", sum(df$warning),
    " skip=", sum(df$skipped), "\n", sep = "")
for (i in seq_along(res)) {
  r <- res[[i]]
  for (x in r$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("--- ", r$test, " [", class(x)[1L], "]\n", sep = "")
      cat(conditionMessage(x), "\n")
    }
  }
}
