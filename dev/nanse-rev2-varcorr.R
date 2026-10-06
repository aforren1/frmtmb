# Reviewer, punch round 1: VarCorr() on a fit with a flat sd block. Its
# delta method multiplies the shown covariance, whose lost rows are NaN,
# so quantities the data determine (another group's sd, the residual
# sd) may come out NaN. y ~ x + (1 | g1) + (1 + x | g2): g1 has a real
# variance, the g2 block is at 0.
#   Rscript dev/nanse-rev2-varcorr.R lane|base|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
# the first seed whose g2 block is flat on the lane build (the same data
# and fit on every build)
for (s in 1:40) {
  set.seed(s)
  n <- 60
  d <- data.frame(x = rnorm(n), g1 = factor(rep(1:12, 5)), g2 = gl(6, 10))
  d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g1] + rnorm(n)
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g1) + (1 + x | g2)),
                              data = d))
  H <- stats::optimHess(fit$opt$par, fit$obj$fn, fit$obj$gr)
  rm <- apply(abs(H), 1, max)
  if (sum(rm < 1e-6) >= 2) break
}
cat("seed", s, "| Hessian rows below 1e-6:",
    paste(ns$outer_par_names(fit)[rm < 1e-6], collapse = " "), "\n")
cat("lost:", paste(names(ns$sdr_of(fit)$se_lost), collapse = " "), "\n")
w <- character()
vc <- withCallingHandlers(VarCorr(fit), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
})
cat(sprintf("VarCorr g1 sd %.4g Est.Error %s | residual sd %.4g Est.Error %s | warnings %d\n",
            vc$g1$sd[1, "Estimate"], format(vc$g1$sd[1, "Est.Error"]),
            vc$residual__$sd[1, "Estimate"],
            format(vc$residual__$sd[1, "Est.Error"]), length(w)))
h <- suppressWarnings(hypothesis(fit, "sd_g1__Intercept > 0", class = NULL))
cat("hypothesis(sd_g1__Intercept) Est.Error", h$hypothesis$Est.Error, "\n")
if (requireNamespace("lme4", quietly = TRUE)) {
  m <- suppressMessages(lme4::lmer(y ~ x + (1 | g1) + (1 + x | g2),
                                   data = d, REML = FALSE))
  print(lme4::VarCorr(m))
}
