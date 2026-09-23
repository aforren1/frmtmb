# Lane wt-correct, item 7: frmtmb's REML = TRUE against mgcv's REML on a
# location-scale smooth. frmtmb integrates the mu coefficients only;
# mgcv's LAML integrates every coefficient of both predictors. This
# measures the gap, and the ML gap on the same data as the control.
# Usage: Rscript dev/correct-reml-mgcv.R <arm>   (base | lane)
arm <- commandArgs(trailingOnly = TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "mgcv",
    format(packageVersion("mgcv")), "\n")

one <- function(seed, n = 500) {
  set.seed(seed)
  dd <- data.frame(x = runif(n), z = runif(n))
  dd$y <- rnorm(n, sin(3 * dd$x), exp(0.8 * cos(5 * dd$z) - 0.5))
  out <- list()
  for (m in c("ML", "REML")) {
    fit <- frm(bf(y ~ s(x), sigma ~ s(z)) + gaussian(), data = dd,
               REML = identical(m, "REML"))
    ref <- mgcv::gam(list(y ~ s(x), ~ s(z)), data = dd,
                     family = mgcv::gaulss(b = 0), method = m)
    mu_f <- as.numeric(frm_linpred(fit, type = "response"))
    sg_f <- as.numeric(frm_linpred(fit, dpar = "sigma", type = "response"))
    mu_g <- fitted(ref)[, 1]
    sg_g <- 1 / fitted(ref)[, 2]
    # smoothing parameters: mgcv's sp against 1 / sd^2 of frmtmb's
    # penalized block, the same penalty scale (gaulss has scale 1)
    sds <- unlist(lapply(fit$frame$re_blocks, function(b) {
      exp(fit$estimates$theta[b$theta_idx])
    }))
    lam_f <- 1 / sds^2
    out[[m]] <- c(
      mu_gap = max(abs(mu_f - mu_g)) / diff(range(mu_g)),
      sg_gap = max(abs(sg_f - sg_g)) / diff(range(sg_g)),
      logsp_x = log(lam_f[1]) - log(ref$sp[1]),
      logsp_z = log(lam_f[2]) - log(ref$sp[2]))
  }
  out
}
res <- lapply(c(44, 45, 46, 47, 48), function(s) {
  r <- one(s)
  cat(sprintf("seed %d  ML:   mu %.4g sigma %.4g logsp x %.4g z %.4g\n",
              s, r$ML[1], r$ML[2], r$ML[3], r$ML[4]))
  cat(sprintf("seed %d  REML: mu %.4g sigma %.4g logsp x %.4g z %.4g\n",
              s, r$REML[1], r$REML[2], r$REML[3], r$REML[4]))
  r
})
