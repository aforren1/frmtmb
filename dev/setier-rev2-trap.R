# Reviewer of lane setier, re-check: fits that stop at a group sd near
# zero although lme4's maximum has it well inside (dev/setier-rev2-scale2.R:
# ri20 seeds 19 and 36 with y x 1e-3, seeds 21 and 23 with y x 1e3).
# What frmtmb tells the user, and the log-likelihood it leaves.
#   Rscript dev/setier-rev2-trap.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
for (cs in list(c(19, 1e-3), c(36, 1e-3), c(21, 1e3), c(23, 1e3))) {
  set.seed(cs[1])
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * cs[2]
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  cat(sprintf("seed %g y x %g: code %d (%s) | logLik - lme4 %.4f | sd/sigma %.3g (lme4 %.3g)\n",
              cs[1], cs[2], f$opt$convergence, f$opt$message,
              as.numeric(logLik(f)) - as.numeric(logLik(l4)),
              exp(f$estimates$theta[1]) / exp(f$estimates$betad[1]),
              attr(VarCorr(l4)$g, "stddev") / sigma(l4)))
  for (x in w) cat("   W:", substr(x, 1, 110), "\n")
  for (x in m) cat("   M:", substr(x, 1, 110), "\n")
}
