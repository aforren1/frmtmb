# Lane wt-mvprior: one test file, one R process. dev/release/run-tests.R
# with the lane library first, plus an optional trace.
#
# With MVPRIOR_TRACE=<file>, every time one of this lane's refusals fires (a
# multivariate prior without resp, a location prior without nlpar or dpar,
# class b with no slope, resp where none applies, and the sd refusal a
# several-location family now meets) one line is appended to <file>: the test
# file and the message. A refusal a test swallows with try() or expect_error()
# still logs, so the log lists every call site the suite reaches, not only the
# ones that fail.
#
# The trace is an instrument, not a tier: it wraps frm_stop(), which
# changes the call that test-conditions.R reads back and fails one
# block of test-importance.R (both pass without it, on both arms:
# dev/mvprior-log/detect-suite.log against the reruns in the findings).
# Tier counts come from runs WITHOUT MVPRIOR_TRACE.
#
# Usage: Rscript dev/mvprior-run-tests.R <package> <path-to-test-file>
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))

suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
suppressMessages(library(p, character.only = TRUE))

tr <- Sys.getenv("MVPRIOR_TRACE")
if (nzchar(tr)) {
  assign("mvprior_trace_file", tr, envir = globalenv())
  assign("mvprior_test_file", basename(f), envir = globalenv())
  invisible(suppressMessages(trace(
    "frm_stop", where = asNamespace("frmtmb"), print = FALSE,
    tracer = quote({
      .m <- tryCatch(paste0(...), error = function(e) "")
      if (grepl(paste0("names no parameter of a multivariate|",
                       "location is nonlinear|several distributional|",
                       "no population-level slope|only to a multivariate|",
                       "takes no resp|No random-effect SDs"), .m)) {
        cat(get("mvprior_test_file", envir = globalenv()), "\t",
            gsub("[[:space:]]+", " ", substr(.m, 1, 160)), "\n",
            file = get("mvprior_trace_file", envir = globalenv()),
            append = TRUE, sep = "")
      }
    }))))
}

r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })

if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      "\n", sep = "")
  bad <- r[r$failed > 0 | r$error, c("test"), drop = FALSE]
  if (nrow(bad)) {
    cat(paste0("FAILED ", basename(f), ": ", bad$test, "\n"), sep = "")
  }
}
