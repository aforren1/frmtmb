# Lane ceplot: one test file, one R process, with the package attached.
#   Rscript dev/ceplot-runtest.R <lane|base> <package> <test file>
# lane: the lane's library first, then rellib-r4 (the other extensions
# of 0.66.0), then the user library. base: rellib-r4 alone, the 0.66.0
# build every new test is seen failing on.
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
p <- a[2]
f <- a[3]
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package(p), "| frmtmb:", find.package("frmtmb"), "\n")
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", arm, " ", p, " ", basename(f), " LOADERROR ",
        conditionMessage(e), "\n", sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", arm, " ", p, " ", basename(f), ": tests=", sum(r$nb),
      " failed=", sum(r$failed), " error=", sum(r$error), " skipped=",
      sum(r$skipped), " warning=", sum(r$warning), " passed=",
      sum(r$passed), "\n", sep = "")
  bad <- r$test[r$failed > 0 | r$error]
  for (t in bad) cat("  fails:", t, "\n")
  w <- r$test[r$warning > 0]
  for (t in w) cat("  warns:", t, "\n")
}
