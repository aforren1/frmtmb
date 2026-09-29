## One test file, one R process. FRMTMB_LIB=base runs it against the
## 0.64.0 reference build, which is how a pin is SEEN TO FAIL.
## Usage: Rscript dev/thresrefit-runtest.R <test-file.R>
arm <- Sys.getenv("FRMTMB_LIB", "lane")
lib <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
args <- commandArgs(trailingOnly = TRUE)
f <- args[1L]
root <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit"
path <- file.path(root, "tests/testthat", f)
library(frmtmb)
library(testthat)
setwd(file.path(root, "tests/testthat"))
res <- test_file(path, package = "frmtmb",
                 env = test_env("frmtmb"), reporter = "silent")
df <- as.data.frame(res)
cat("RESULT", f, "arm=", arm,
    " pass=", sum(df$passed),
    " fail=", sum(df$failed),
    " err=", sum(df$error),
    " skip=", sum(df$skipped), "\n", sep = "")
## the failures themselves, so a "seen to fail" record has the message
for (i in seq_along(res)) {
  r <- res[[i]]
  for (x in r$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("--- ", r$test, " [", class(x)[1L], "]\n", sep = "")
      cat(conditionMessage(x), "\n")
    }
  }
}
