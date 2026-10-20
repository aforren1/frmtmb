# Lane adhess: frmtmb.coupling test-coherence.R at eta = 16, the loop
# step that fails on the lane. What the fit's covariance and the
# coherence interval are on each arm, against the exact AD Hessian.
#   ADHESS_ARM=lane|base Rscript dev/adhess-coh.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
suppressMessages(library(frmtmb.coupling))
tf <- "C:/Users/adf44/source/r/frmtmb-wt-adhess/extensions/frmtmb.coupling/tests/testthat/test-coherence.R"
ex <- parse(tf)
for (e in ex[1:3]) eval(e)
set.seed(201)
for (eta in c(-3, 0, 3, 16)) {
  d <- cp_units(rep(eta, 15L), rep(0.4, 15L), n = 4L)
  r <- cap(frm(cp_bf(), family = cross_wishart(), data = d))
  fit <- r$v
  co <- try(frm_coherence(fit, newdata = d[1, ]), silent = TRUE)
  p <- fit$opt$par
  H <- fit$obj$he(p)
  ev <- eigen((H + t(H)) / 2, only.values = TRUE)$values
  Vf <- suppressWarnings(vcov(fit, full = TRUE))
  cat(sprintf("eta %g | code %d | par %s | eig(AD H) %s | warn %s | lost %s\n",
              eta, fit$opt$convergence, paste(signif(p, 5), collapse = ","),
              paste(signif(ev, 3), collapse = ","), wtag(r$w), lost_of(fit)))
  cat("  SE:", format(sqrt(diag(Vf)), digits = 5), "\n")
  if (inherits(co, "try-error")) cat("  coherence: error",
                                    conditionMessage(attr(co, "condition")), "\n")
  else cat(sprintf("  coherence %.12f [%.12f, %.12f] eta %.6f\n",
                   co$.estimate, co$.lower, co$.upper, co$.eta))
}
