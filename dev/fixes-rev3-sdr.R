# Reviewer of lane fixes, final check: on a smooth fit, whose par_units
# now come from smooth_fx_units(), autoscale_sdreport()'s transformed
# covariance against RTMB::sdreport()'s plain one.
#   Rscript dev/fixes-rev3-sdr.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (seed in 1:3) {
  set.seed(seed)
  d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
  fit <- suppressWarnings(frm(bf(y ~ s(x1) + s(x2) + x0), data = d))
  V <- unname(vcov(fit, full = TRUE))
  P <- unname(RTMB::sdreport(fit$obj)$cov.fixed)
  cat(sprintf("seed %d units %s: max |V - P| / max |P| = %.3g\n", seed,
              paste(format(fit$par_units, digits = 3), collapse = ","),
              max(abs(V - P)) / max(abs(P))))
}
