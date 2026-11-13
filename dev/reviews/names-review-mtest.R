# Review of lane names: one frmtmb.sample test file on the names +
# sampfront merged build (dev/reviews/names-review-mtree.sh).
# Usage: Rscript dev/reviews/names-review-mtest.R <file> [gated]
a <- commandArgs(TRUE)
SP <- "C:/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad"
.libPaths(c(file.path(SP, "mlib"), "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if ("gated" %in% a) Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
p <- "frmtmb.sample"
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package("frmtmb"), find.package(p), "\n")
f <- file.path(SP, "mtree", p, "tests/testthat", a[1])
tr <- testthat::test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent", stop_on_failure = FALSE)
for (t in tr) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) {
      cat("FAILURE [", t$test, "]: ", substr(conditionMessage(e), 1, 600),
          "\n", sep = "")
    }
  }
}
r <- as.data.frame(tr)
cat(sprintf("RESULT %s pass=%d fail=%d err=%d skip=%d\n", a[1],
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped)))
