# Reviewer of lane fixes, final check: the ported update of
# brmsfit-methods:955 (fixture without priors) from a usable start.
#   Rscript dev/fixes-rev3-u955.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit2 <- brms_fixture(2)
d <- brms_fixture_data(2)
r <- tryCatch(withCallingHandlers(
  update(fit2, formula. = bf(count ~ a + b, nl = TRUE),
         start = list(beta = c(mean(d$count) / 2, 0, mean(d$count) / 2, 0))),
  warning = function(w) {
    cat("[warning]", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) {cat("ERROR", conditionMessage(e), "\n"); NULL})
if (!is.null(r)) {
  cat("conv", r$opt$convergence, " logLik", format(as.numeric(logLik(r)),
                                                    digits = 10), "\n")
  print(suppressWarnings(fixef(r))[, 1:2])
}
