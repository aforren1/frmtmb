# Reviewer runner: one test file, one R process, package attached.
# Usage: Rscript formrobust-rev-run-tests.R <package> <test file>
# REVLIB="" gives the base arm (rellib-r4 only).
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
Sys.setenv(R_MAKEVARS_USER = mk)
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(testthat))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]; f <- a[2]
mut <- Sys.getenv("REVMUTANT", "")
suppressMessages(library(p, character.only = TRUE))
if (nzchar(mut)) {
  # a mutant: a script that reassigns functions in the namespace
  source(mut, local = new.env())
  cat("MUTANT ", mut, "\n", sep = "")
}
res <- tryCatch(
  test_file(f, package = p, env = testthat::test_env(p),
            reporter = "silent"),
  error = function(e) {
    cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
        sep = "")
    NULL
  })
r <- if (!is.null(res)) as.data.frame(res)
if (!is.null(res)) {
  for (tt in res) {
    for (ex in tt$results) {
      cl <- class(ex)[1]
      if (cl %in% c("expectation_failure", "expectation_error",
                    "expectation_warning")) {
        cat("DETAIL [", cl, "] ", tt$test, ": ",
            substr(conditionMessage(ex), 1, 400), "\n", sep = "")
      }
    }
  }
}
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), " tests=", nrow(r), "\n", sep = "")
}
cat("lib: ", find.package(p), " frmtmb=", find.package("frmtmb"), "\n",
    sep = "")
