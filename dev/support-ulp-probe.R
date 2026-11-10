# Lane support: do the new trunc(lb <= 0) assertions of test-cens-trunc.R
# hold on every rounding? The test's data (seed 401) scaled by
# (1 + 2e-15 k), k = 1..20, and for each family and lb in {-1, 0} the
# two quantities the test bounds: |logLik difference| / (64 eps |logLik|)
# and max |estimate difference| / (1e-6 max SE). Both must stay <= 1.
#
#   Rscript dev/support-ulp-probe.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", paste0(
  "C:/Users/adf44/source/r/wt-support-lib;",
  "C:/Users/adf44/source/r/rellib-r19"))
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(401)
e <- stats::rexp(300)
ys <- list(exponential = 2 * e, weibull = 2 * e,
           frechet = 2 / gamma(1 - 1 / 3) * e^(-1 / 3))
worst <- c(ll = 0, est = 0)
for (k in 1:20) {
  for (fam in list(exponential(), weibull(), frechet())) {
    d <- data.frame(y = ys[[fam$family]] * (1 + 2e-15 * k))
    f0 <- frm(y ~ 1, family = fam, data = d)
    for (lb in c(-1, 0)) {
      d$lb <- lb
      f1 <- frm(y | trunc(lb = lb) ~ 1, family = fam, data = d)
      l0 <- as.numeric(logLik(f0))
      r1 <- abs(as.numeric(logLik(f1)) - l0) /
        (64 * .Machine$double.eps * abs(l0))
      r2 <- max(abs(fixef(f1)[, 1] - fixef(f0)[, 1])) /
        (1e-6 * max(fixef(f0)[, 2]))
      worst <- pmax(worst, c(r1, r2))
    }
  }
}
cat(sprintf(paste("draws 20 x 3 families x 2 bounds: worst logLik ratio",
                  "%.3g, worst estimate ratio %.3g\n"), worst[1], worst[2]))
