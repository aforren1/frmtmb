# The Jacobian mutants send every fit up the refusal's ladder, so the whole file
# takes very long. Run named tests on them, one at a time.
# Usage: Rscript dev/review-impsamp-jac.R <mutant>
m <- commandArgs(TRUE)[[1]]
LIB <- "C:/Users/adf44/source/r/rev-impsamp-lib"
.libPaths(c(file.path(LIB, m), "C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
stopifnot(startsWith(find.package("frmtmb"), LIB))
f <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/tests/testthat/test-importance.R"
for (d in c("the objective's curvature along an intercept is not a sample",
            "a gaussian fit stays at the exact optimum at any draw count",
            "gaussian: the corrected objective at its anchor IS Laplace's",
            "the corrected objective's gradient matches numDeriv",
            "the textbook binomial settles at the draw counts users pick",
            "scalar intercept agrees with quadrature and glmer(nAGQ = 25)")) {
  r <- as.data.frame(test_file(f, package = "frmtmb", env = test_env("frmtmb"),
                               reporter = "silent", desc = d, stop_on_failure = FALSE))
  cat(sprintf("%s | %s: pass=%d fail=%d error=%d\n", m, d, sum(r$passed), sum(r$failed),
              sum(r$error)))
}
