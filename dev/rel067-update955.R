# brmsfit-methods:955 on the release build and on the lane builds:
# update(fit2, formula. = bf(count ~ a + b, nl = TRUE)) on the ported
# suite's fixture 2 (brms's brmsfit_example2).
#
#   Rscript dev/rel067-update955.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", lib, as.character(packageVersion("frmtmb")), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit2 <- brms_fixture(2)
cat("fit2 formula:\n"); print(formula(fit2))
up <- tryCatch(withCallingHandlers(
  update(fit2, formula. = bf(count ~ a + b, nl = TRUE)),
  message = function(m) {
    cat("  [message]", conditionMessage(m)); invokeRestart("muffleMessage")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  })
if (!is.null(up)) {
  cat("class:", class(up), "\n"); print(formula(up))
  cat("logLik:", format(as.numeric(logLik(up)), digits = 10), "\n")
}
