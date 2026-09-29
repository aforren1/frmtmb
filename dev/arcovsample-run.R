# Lane wt-arcovsample: one test file, one R process.
#
#   Rscript dev/arcovsample-run.R <lib> <package> <test-file>
#
# <lib> is "lane" (this lane's private library) or "ref" (the read-only
# reference build of the base commit, for a pin seen to fail).

a <- commandArgs(trailingOnly = TRUE)
LIB <- switch(a[1L],
              lane = "C:/Users/adf44/source/r/wt-arcovsample-lib",
              ref = "C:/Users/adf44/source/r/rellib-r3",
              stop("first argument must be 'lane' or 'ref'"))
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))

suppressMessages(library(testthat))
p <- a[2L]
f <- a[3L]
suppressMessages(library(p, character.only = TRUE))

r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", a[1L], " ", p, " ", basename(f), " LOADERROR ",
        conditionMessage(e), "\n", sep = "")
    NULL
  })

if (!is.null(r)) {
  cat("RESULT ", a[1L], " ", p, " ", basename(f), " pass=", sum(r$passed),
      " fail=", sum(r$failed), " err=", sum(r$error), " skip=",
      sum(r$skipped), "\n", sep = "")
  bad <- r[r$failed > 0 | r$error > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat("  BAD [", bad$test[i], "] fail=", bad$failed[i], " err=",
        bad$error[i], "\n", sep = "")
  }
}
