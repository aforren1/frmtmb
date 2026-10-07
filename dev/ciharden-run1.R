# One test file, one R process, for the fragility scan
# (dev/ciharden-scan.sh). The same runner as dev/release/run-tests.R,
# with two differences the scan needs: the library stack is an argument,
# so a lane build and the base build run through one script; and every
# non-passing expectation is printed with its source line, so that two
# configurations can be compared expectation by expectation rather than
# by counts that can cancel.
#
# Usage: Rscript dev/ciharden-run1.R <lib dir or "base"> <package> <file>
local({
# ROUND_BASE lets a later round rerun the scan against its own base
base <- Sys.getenv("ROUND_BASE", "C:/Users/adf44/source/r/rellib-r7")
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1], base, user))
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
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
# the BLAS actually loaded: a 600 x 600 product is about 0.1 s with the
# reference BLAS and under 0.01 s with OpenBLAS
cat("blas: ", R.home(), " threads=", Sys.getenv("OPENBLAS_NUM_THREADS"),
    " gemm600=", system.time({
      m <- matrix(1, 600, 600); m %*% m
    })[["elapsed"]], "\n", sep = "")
res <- NULL
r <- tryCatch({
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})
if (!is.null(res)) {
  for (t in res) {
    for (x in t$results) {
      k <- if (inherits(x, "expectation_failure")) "FAIL" else
        if (inherits(x, "expectation_error")) "ERROR" else
          if (inherits(x, "expectation_skip")) "SKIP" else
            if (inherits(x, "expectation_warning")) "WARN" else NA
      if (!is.na(k)) {
        ln <- if (!is.null(x$srcref)) {
          paste0(basename(utils::getSrcFilename(x$srcref)), ":",
                 x$srcref[1])
        } else "?"
        cat("DETAIL\t", k, "\t", ln, "\t", t$test, "\t",
            substr(gsub("[\r\n\t]+", " | ", conditionMessage(x)), 1, 400),
            "\n", sep = "")
      }
    }
    # a test's pass count can move with the platform when it loops over
    # what a fit returns, so the scan compares it too
    cat("TEST\t", t$test, "\t",
        sum(vapply(t$results, inherits, logical(1),
                   "expectation_success")), "\n", sep = "")
  }
}
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})
