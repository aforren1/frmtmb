# Punch round 2, RB2 fixture: where does summary()'s "NaNs produced"
# come from, and what do ranef(condVar) and conditional_effects(re_formula
# = NULL) give on the singular g2 block?
.libPaths(c(if (!identical(commandArgs(TRUE)[1], "base")) "C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), g1 = factor(rep(1:12, 5)), g2 = gl(6, 10))
d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g1] + rnorm(n)
fit <- suppressWarnings(frm(bf(y ~ x + (1 | g1) + (1 + x | g2)), data = d))
sdr <- frmtmb:::sdr_of(fit)
cat("diag.cov.random range:", range(sdr$diag.cov.random), " negative:",
    sum(sdr$diag.cov.random < 0), "of", length(sdr$diag.cov.random), "\n")
options(warn = 2)
r <- tryCatch(summary(fit), error = function(e) {
  cat("summary error/warning:", conditionMessage(e), "\n")
  print(sys.calls())
  NULL
})
tb <- tryCatch(withCallingHandlers(summary(fit), warning = function(w) {
  cat("WARNING:", conditionMessage(w), "\n")
  calls <- sys.calls()
  for (cl in utils::tail(calls, 12)) cat("  ", deparse(cl)[1], "\n")
  invokeRestart("muffleWarning")
}), error = function(e) NULL)
options(warn = 1)
re <- withCallingHandlers(ranef(fit, condVar = TRUE), warning = function(w) {
  cat("ranef warning:", conditionMessage(w), "\n")
  invokeRestart("muffleWarning")
})
cat("ranef condSD g1:", format(as.numeric(attr(re$g1, "condSD")),
                               digits = 3), "\n")
cat("ranef condSD g2:", format(as.numeric(attr(re$g2, "condSD")),
                               digits = 3), "\n")
ce <- withCallingHandlers(conditional_effects(fit, "x", re_formula = NULL),
  warning = function(w) {
    cat("ce warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
cat("ce re_formula = NULL se__ finite", sum(is.finite(ce[[1]]$se__)), "of",
    nrow(ce[[1]]), "\n")
