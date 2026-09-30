# Reviewer runner (lane ceplot): one test file, one R process, package
# attached. Optional mutant injected into the frmtmb namespace first.
#   Rscript dev/ceplot-rev-runtest.R <lane|base> <package> <file> [mutant]
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
p <- a[2]
f <- a[3]
mut <- if (length(a) >= 4) a[4] else ""
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package(p), "| frmtmb:", find.package("frmtmb"),
    "| cache:", Sys.getenv("FRMTMB_STAN_CACHE"), "\n")
if (nzchar(mut)) {
  source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-mutants2.R")
  apply_mutant(mut)
  cat("mutant:", mut, "\n")
}
r <- tryCatch(
  as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                          reporter = "silent")),
  error = function(e) {
    cat("RESULT ", arm, " ", p, " ", basename(f), " ", mut, " LOADERROR ",
        conditionMessage(e), "\n", sep = "")
    NULL
  })
if (!is.null(r)) {
  cat("RESULT ", arm, " ", p, " ", basename(f), " ", mut, ": tests=",
      sum(r$nb), " failed=", sum(r$failed), " error=", sum(r$error),
      " skipped=", sum(r$skipped), " warning=", sum(r$warning), " passed=",
      sum(r$passed), "\n", sep = "")
  bad <- r$test[r$failed > 0 | r$error]
  for (t in bad) cat("  fails:", t, "\n")
  w <- r$test[r$warning > 0]
  for (t in w) cat("  warns:", t, "\n")
}
