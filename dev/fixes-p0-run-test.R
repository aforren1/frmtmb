# Lane fixes, punch round: dev/fixes-run-test.R with the pre-punch guard.
# A copy of dev/release/run-tests.R that takes the library first and
# prints every failure and error message.
#
#   Rscript dev/fixes-run-test.R <lib> <package> <test file>
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
.libPaths(unique(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(testthat))
p <- a[2]
f <- a[3]
suppressMessages(library(p, character.only = TRUE))
# the pre-punch guard in place of the new one
ns <- asNamespace("frmtmb")
sys.source("C:/Users/adf44/source/r/frmtmb-wt-fixes/dev/fixes-p0-nlcheck.R", envir = ns2 <- new.env(parent = ns))
f0 <- ns2$p0_check_nl_sum_identified
environment(f0) <- ns
assignInNamespace("check_nl_identified", f0, ns = "frmtmb")
cat("guard: pre-punch check_nl_sum_identified
")
cat("lib:", find.package(p), "core:", find.package("frmtmb"), "\n")
res <- tryCatch(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent"),
                error = function(e) {
                  cat("RESULT ", basename(f), " LOADERROR ",
                      conditionMessage(e), "\n", sep = "")
                  NULL
                })
if (!is.null(res)) {
  for (t in res) {
    for (e in t$results) {
      if (inherits(e, c("expectation_failure", "expectation_error",
                        "expectation_warning"))) {
        cat("---", class(e)[1], "in:", t$test, "\n",
            conditionMessage(e), "\n")
      }
    }
  }
  r <- as.data.frame(res)
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
