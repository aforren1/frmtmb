# Reviewer: pre-existing? Beta() on rellib-r3 at high precision (seed 4500).
.libPaths(c("C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (phi in c(500, 2000, 5000, 2e4)) {
  set.seed(4500); y <- rbeta(1000, 0.3 * phi, 0.7 * phi)
  f <- tryCatch(frm(y ~ 1, family = Beta(), data = data.frame(y = y)), error = function(e) e)
  cat(sprintf("phi %g: %s\n", phi, if (inherits(f, "error")) substr(conditionMessage(f), 1, 90) else
    sprintf("ok, phi-hat %.1f", frmtmb:::eval_dpars(f)[[1]]$phi[1])))
}
x <- RTMB::MakeTape(function(p) RTMB::dbeta(0.3, p[1], p[2], log = TRUE), c(1, 1))
for (s in c(1000, 1100, 1500, 2000)) cat(s, x$jacobian(c(0.3 * s, 0.7 * s)), "\n")
for (s in c(1000, 1100, 1500, 2000)) cat("lgamma grad", s, RTMB::MakeTape(function(p) lgamma(p), 1)$jacobian(s), "\n")
