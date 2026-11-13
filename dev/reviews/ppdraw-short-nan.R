# Short check: does predict() leak "NaNs produced" from draw_edge_csd()
# on the cliff family's fits, and which Hessian diagonal is negative?
# Usage: Rscript ppdraw-short-nan.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(5)
dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
a <- stats::rnorm(20, 0, 0.7)[dd$g]
b <- stats::rnorm(20, 0, 0.05)[dd$g] * dd$x
e <- stats::rnorm(200)
dd$y <- 1 + 0.5 * dd$x + a + b + e
fit <- frm(bf(y ~ x + (1 + x | g)), data = dd)
nd <- data.frame(x = 2, g = factor("new"))
set.seed(1)
w <- character(0)
p <- withCallingHandlers(
  predict(fit, newdata = nd, ndraws = 100, allow_new_levels = TRUE),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("M:", substr(conditionMessage(m), 1, 100), "\n")
    invokeRestart("muffleMessage")
  })
cat("warnings from predict():", if (length(w)) w else "(none)", "\n")
hc <- fit$cache$hessian_fixed
if (!is.null(hc$H)) {
  cat("diag(H):\n")
  print(signif(stats::setNames(diag(hc$H), names(fit$opt$par)), 4))
  cat("names in hessian_fixed:", paste(names(hc), collapse = ", "), "\n")
}
