## Reviewer, the risk the worker did not test: helper-fuzz.R MUTES the
## metamorphic invariants on any fit that raises a convergence warning
## (FUZZ_NONCONVERGENCE matches "Large maximum absolute gradient"), so
## removing 40 percent of those warnings puts MORE fits under the
## invariants than before. If any of them fails an invariant, the fuzz
## tier goes red on the lane build and was green on base.
##
## Run with the tier enabled and a bounded plan, one arm per process.
## usage: Rscript gradcheck-rev-17-fuzz.R <core-lib> <size>
a <- commandArgs(TRUE)
CORE <- a[1]
SIZE <- if (length(a) >= 2) a[2] else "120"
.libPaths(c(CORE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true", FRMTMB_FUZZ_N = SIZE)
cat("CORE:", CORE, "  size:", SIZE, "\n")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
FILE <- file.path("C:/Users/adf44/source/r/frmtmb-wt-gradcheck",
                  "tests/testthat/test-fuzz.R")
tenv <- new.env(parent = asNamespace("frmtmb"))
res <- try(testthat::test_file(FILE, reporter = "silent", env = tenv),
           silent = TRUE)
if (inherits(res, "try-error")) {
  cat("RESULT test-fuzz.R ABORTED\n"); cat(as.character(res), "\n")
} else {
  df <- as.data.frame(res)
  cat(sprintf("RESULT test-fuzz.R pass=%d fail=%d error=%d skip=%d\n",
              sum(df$passed), sum(df$failed), sum(df$error),
              sum(df$skipped)))
  for (tb in res) for (ex in tb$results) {
    if (inherits(ex, "expectation_failure") ||
        inherits(ex, "expectation_error") ||
        inherits(ex, "expectation_skip")) {
      cat("### ", class(ex)[1], ": ", tb$test, "\n", sep = "")
      cat(conditionMessage(ex), "\n")
    }
  }
}
