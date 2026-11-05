# Ledger row brmsfit-methods:955: update(fit2, bf(count ~ a + b, nl = TRUE))
# stopped at its zero start on Gamma("identity"), which nl_domain_start() now
# moves. The row's own call, on the suite's fixture.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit2 <- brms_fixture(2)
up <- tryCatch(withCallingHandlers(
  update(fit2, formula. = bf(count ~ a + b, nl = TRUE)),
  message = function(m) cat("message:", substr(conditionMessage(m), 1, 160), "\n"),
  warning = function(w) {
    cat("warning:", substr(conditionMessage(w), 1, 160), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) e)
if (inherits(up, "error")) {
  cat("ERROR:", substr(conditionMessage(up), 1, 300), "\n")
} else {
  cat("logLik", as.numeric(logLik(up)), "code", up$opt$convergence, "\n")
  print(round(fixef(up)[, 1], 4))
}
