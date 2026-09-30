# Reviewer of lane defects: one test file, one process, with the gates
# open, on the lane build (arm "lane") or the base build (arm "base").
# Prints a RESULT line and every failure, error and escaped warning with
# its location.
#   Rscript dev/defects-rev-runtest.R <arm> <pkg> <file>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]; p <- a[2]; f <- a[3]
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p,
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
res <- test_file(f, package = p, env = test_env(p), reporter = "list",
                 stop_on_failure = FALSE)
r <- as.data.frame(res)
cat("RESULT", arm, p, format(packageVersion(p)),
    dirname(getNamespaceInfo(p, "path")), basename(f),
    "pass=", sum(r$passed), "fail=", sum(r$failed), "err=",
    sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning), "\n")
for (t in res) {
  for (e in t$results) {
    if (inherits(e, c("expectation_failure", "expectation_error",
                      "expectation_warning", "expectation_skip"))) {
      sr <- e$srcref
      loc <- if (is.null(sr)) "?" else as.character(sr[1])
      cat("--", class(e)[1], "line", loc, "in:", t$test, "\n  ",
          gsub("\n", "\n   ", substr(conditionMessage(e), 1, 500)), "\n")
    }
  }
}
