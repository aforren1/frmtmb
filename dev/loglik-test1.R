# Lane loglik: one test file, one process, against the lane's library.
#   Rscript dev/loglik-test1.R <package> <test file> [base]
# "base" drops the lane's library (the "before" arm on rellib-r13).
a <- commandArgs(trailingOnly = TRUE)
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (length(a) >= 3L && a[3] == "base") libs <- libs[-1]
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE =
               "C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/stan-cache")
}
p <- a[1]
f <- a[2]
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), p, dirname(find.package(p)), "\n")
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "summary", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), "\n", sep = "")
