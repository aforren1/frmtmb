# One test file, one R process, for lane ordinal: dev/release/run-tests.R
# with the lane's private library first (arm "lane") or the round's
# base build alone (arm "base"). The base library stays behind the lane
# library, so an extension this lane did not change loads the lane's
# frmtmb.
#
# Usage: Rscript dev/ordinal-run-tests.R <arm> <package> <test file>

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
p <- a[2]
f <- a[3]
LIB <- c("C:/Users/adf44/source/r/rellib-r4",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) LIB <- c("C:/Users/adf44/source/r/wt-ordinal-lib", LIB)
.libPaths(LIB)
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package(p), " frmtmb:", find.package("frmtmb"), "\n")
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
      " warn=", sum(r$warning), "\n", sep = "")
  bad <- r[r$failed > 0 | r$error, c("test", "failed", "error"),
           drop = FALSE]
  if (nrow(bad)) {
    cat("FAILING TESTS:\n")
    print(bad, row.names = FALSE)
  }
}
