# Lane wt-arcovsample: one test file, one R process, with the brms fit
# tier ON. Stan compiles fresh on this machine (no cache directory
# exists), so this is minutes per model.
#
#   Rscript dev/arcovsample-run-gated.R <package> <test-file>
a <- commandArgs(trailingOnly = TRUE)
# ARCOVSAMPLE_REF=true measures the base commit instead, for a pin seen
# to fail. That library is READ-ONLY and nothing here writes to it.
LIB <- if (identical(Sys.getenv("ARCOVSAMPLE_REF"), "true")) {
  "C:/Users/adf44/source/r/rellib-r3"
} else {
  "C:/Users/adf44/source/r/wt-arcovsample-lib"
}
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))

suppressMessages(library(testthat))
p <- a[1L]
f <- a[2L]
suppressMessages(library(p, character.only = TRUE))
r <- test_file(f, package = p, env = testthat::test_env(p),
               reporter = "silent")
d <- as.data.frame(r)
cat("RESULT gated ", p, " ", basename(f), " pass=", sum(d$passed),
    " fail=", sum(d$failed), " err=", sum(d$error), " skip=",
    sum(d$skipped), "\n", sep = "")
for (tr in r) {
  for (x in tr$results) {
    if (inherits(x, "expectation_skip")) {
      cat("SKIP [", tr$test, "] ", conditionMessage(x), "\n", sep = "")
    }
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("BAD [", tr$test, "] ", conditionMessage(x), "\n", sep = "")
    }
  }
}
