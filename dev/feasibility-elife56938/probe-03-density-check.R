# Probe 3: is wiener()'s density the same object RWiener computes, and how do the
# fitted coefficients come back on the natural scale?
#
# The M1 recovery run returned a drift estimate a few standard errors from the
# simulating value while the boundary separation landed exactly, so the density itself
# needs checking before any recovery claim is made.

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.ddm); library(RWiener)
})

cat("=== density: frmtmb.ddm vs RWiener (ground truth) ===\n")
grid <- expand.grid(dt = c(0.02, 0.05, 0.2, 0.6, 1.5), v = c(0, 1, 3, 5.4),
                    a = c(0.8, 1.6, 3.6), up = c(0L, 1L))
lp_pkg <- frmtmb.ddm:::ddm_lpdf_both(grid$dt, grid$v, grid$a, rep(0.5, nrow(grid)), grid$up)
TAU <- 0.1   # RWiener rejects tau = 0, so shift both the argument and the offset
lp_rw <- mapply(function(dt, v, a, up)
  RWiener::dwiener(dt + TAU, alpha = a, tau = TAU, beta = 0.5,
                   delta = v, resp = if (up == 1L) "upper" else "lower", give_log = TRUE),
  grid$dt, grid$v, grid$a, grid$up)
grid$u <- grid$dt / grid$a^2
err <- lp_pkg - lp_rw
cat("max abs log-density error:", format(max(abs(err)), digits = 4), "\n")
worst <- order(-abs(err))[1:8]
print(cbind(grid[worst, c("dt", "v", "a", "up", "u")],
            pkg = round(lp_pkg[worst], 5), rwiener = round(lp_rw[worst], 5),
            err = signif(err[worst], 3)))
cat("\nerror by u decile (u = t / a^2, the series-blend variable):\n")
print(tapply(abs(err), cut(grid$u, breaks = quantile(grid$u, 0:5 / 5), include.lowest = TRUE),
             max))

cat("\n=== coefficient extraction and natural-scale back-transform ===\n")
set.seed(1)
COH <- c(0, 0.032, 0.064, 0.128, 0.256, 0.512)
d <- do.call(rbind, lapply(COH, function(cc) {
  x <- RWiener::rwiener(440, alpha = 1.6, tau = 0.28, beta = 0.5, delta = 10.5 * cc)
  data.frame(rt = x$q, upper = as.integer(x$resp == "upper"), coh = cc,
             cohf = factor(cc, levels = COH))
}))
m <- frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = d)
cat("class(fixef):", class(fixef(m)), "\n")
print(fixef(m))
cat("\ncoef():\n"); print(coef(m))
cat("\nnatural scale via fitted():\n")
for (dp in c("mu", "bs", "ndt")) cat(" ", dp, "=", signif(fitted(m, dpar = dp)[1], 6), "\n")
cat("min(rt) =", min(d$rt), " (this caps the ndt link)\n")
