## REVIEW claim 5: run frmtmb.sample's prior test files with the WORKER's
## core build underneath and frmtmb.sample from the reference library, and
## again with both from the reference library, so the label change is
## measured against a control. No install is needed: with the lane library
## first, frmtmb resolves there and frmtmb.sample resolves to rellib-r3.
## Usage: FRMTMB_LIB=[lane|base] Rscript dev/thresrefit-rev-13-sample.R <f>
arm <- Sys.getenv("FRMTMB_LIB", "lane")
core <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(core, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
f <- commandArgs(trailingOnly = TRUE)[1L]
ext <- "C:/Users/adf44/source/r/frmtmb/extensions/frmtmb.sample"
library(frmtmb)
library(frmtmb.sample)
library(testthat)
cat("core frmtmb from:",
    dirname(dirname(getNamespaceInfo("frmtmb", "path"))), "\n")
cat("frmtmb.sample from:",
    dirname(dirname(getNamespaceInfo("frmtmb.sample", "path"))), "\n")
setwd(file.path(ext, "tests/testthat"))
res <- test_file(file.path(ext, "tests/testthat", f),
                 package = "frmtmb.sample",
                 env = test_env("frmtmb.sample"), reporter = "silent")
df <- as.data.frame(res)
cat("SAMPLERESULT ", f, " arm=", arm,
    " blocks=", nrow(df), " pass=", sum(df$passed),
    " fail=", sum(df$failed), " err=", sum(df$error),
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
