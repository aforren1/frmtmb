# Reviewer of lane loglik: one test file, one process.
#   Rscript dev/review-loglik-test1.R <package> <file> <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
base <- identical(a[3], "base")
lib <- if (base) character() else "C:/Users/adf44/source/r/wt-loglik-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/stan-cache")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", find.package("frmtmb"), find.package(p), "\n")
r <- as.data.frame(test_file(f, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
w <- sum(vapply(r$result, function(x) {
  sum(vapply(x, inherits, TRUE, "expectation_warning"))
}, 0))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", w, "\n", sep = "")
bad <- r[r$failed > 0 | r$error, "test"]
if (length(bad)) cat("FAILING BLOCK:", bad, sep = "\n  ")
for (x in r$result) for (e in x) {
  if (inherits(e, c("expectation_failure", "expectation_error"))) {
    cat("MSG:", substr(conditionMessage(e), 1, 300), "\n")
  }
}
