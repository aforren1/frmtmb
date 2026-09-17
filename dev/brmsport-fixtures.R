# Build the frmtmb stand-ins for brms's brmsfit_example1..6 and report
# what each costs and whether it converged. The fixtures themselves are
# defined in tests/testthat/helper-brms-suite.R (brms_fixture_spec());
# this script only exercises them.
#
#   Rscript dev/brmsport-fixtures.R > dev/brmsport-log/fixtures.txt 2>&1
lib <- "C:/Users/adf44/source/r/brmsport-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
for (k in 1:6) {
  cat("\n== fixture", k, "\n")
  spec <- brms_fixture_spec(k)
  print(spec$formula)
  cat("brms formula:", spec$brms, "\n")
  cat("changed:", spec$changed, "\n")
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(brms_fixture(k), error = function(e) e)
  t1 <- proc.time()[["elapsed"]]
  if (inherits(fit, "error")) {
    cat("ERROR:", conditionMessage(fit), "\n")
    next
  }
  cat(sprintf("fitted in %.1f s; nobs %d; convergence %s; logLik %.6f\n",
              t1 - t0, nobs(fit), format(fit$opt$convergence),
              as.numeric(logLik(fit))))
  w <- attr(fit, "brmsport_warnings")
  if (length(w)) cat("warnings:", unique(w), sep = "\n  ")
  cat("variables:", head(variables(fit), 30), "\n")
}
