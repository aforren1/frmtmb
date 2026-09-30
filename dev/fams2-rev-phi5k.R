# Reviewer: end to end, xbeta (kappa free) vs Beta() at phi 3e3..1e4
# (seed 4200 construction of dev/fams2-rev-dens.R section B).
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (phi in c(3e3, 5e3, 1e4)) {
  set.seed(4200); n <- 1000; mu <- 0.047; kap <- 0.05
  z <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
  f <- tryCatch(frm(y ~ 1, family = xbeta(), data = data.frame(y = y)), error = function(e) e)
  fb <- tryCatch(frm(y ~ 1, family = Beta(), data = data.frame(y = y[y > 0])), error = function(e) e)
  cat(sprintf("phi %.0e (zeros %d): xbeta %s | Beta() on the interior rows %s\n", phi, sum(y == 0),
    if (inherits(f, "error")) paste("ERROR", substr(conditionMessage(f), 1, 70)) else sprintf("ok phi-hat %.0f", frmtmb:::eval_dpars(f)[[1]]$phi[1]),
    if (inherits(fb, "error")) "ERROR" else sprintf("ok phi-hat %.0f", frmtmb:::eval_dpars(fb)[[1]]$phi[1])))
}
