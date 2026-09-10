# NIT 9: the smoothness of the correction FACTOR itself, measured on
# the shape rather than on an objective built over it.
#
# Round 1's probe swept `k21` through a solved model, and at the
# shipped `n_ss` = 20 a ratio of 0.05 leaves the run-in 1e-26 from its
# limit, so the damping shuts the gate on both sides and the low corner
# was never in the code path. The reviewer's construction removes the
# model: with ya = 0, yb = 1/r and yc = 1/r + 1 the differences are
# 1/r and 1, the damping is negligible, and the ratio the shipped
# function measures IS r, so the correction it returns is F(r) times 1.
#
# The reading is the ORDER of the junction: gap ~ eps^k means C(k).
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-26-corners.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()

RCAP <- frmtmb.ode:::ode_ss_rcap
RLOW <- frmtmb.ode:::ode_ss_rlow
FF <- function(r) {
  # 1 / r is the whole construction, so r = 0 exactly is the one point
  # it cannot express; the correction there is 0 by inspection, since
  # the gate is 0 and the factor carries r.
  if (r == 0) return(c(f = 0, r = 0))
  z <- frmtmb.ode:::ode_ss_extrapolate(0, 1 / r, 1 / r + 1, 1e-300,
                                       1e-300)
  c(f = as.numeric(z[["y"]]) - (1 / r + 1), r = as.numeric(z[["r"]]))
}
cat("\nthe ratio the function reads against the r asked for:\n")
for (r in c(0.02, RLOW, 0.5, RCAP, 1, 1.2))
  cat(sprintf("  asked %-8g read %-16.12g  error %.2e\n", r,
              FF(r)[["r"]], abs(FF(r)[["r"]] - r)))

one_sided <- function(r0, eps) {
  f0 <- FF(r0)[["f"]]
  lo <- (f0 - FF(r0 - eps)[["f"]]) / eps
  hi <- (FF(r0 + eps)[["f"]] - f0) / eps
  c(lo = lo, hi = hi, gap = abs(hi - lo))
}
cat("\n== the four junctions ==\n")
cat(sprintf("%-26s %11s %11s %11s %11s\n", "junction", "eps 1e-2",
            "1e-3", "1e-4", "1e-5"))
for (z in list(list("r = 0", 0), list("r = rlow", RLOW),
               list("r = rcap", RCAP), list("r = 1", 1))) {
  g <- vapply(c(1e-2, 1e-3, 1e-4, 1e-5),
              function(e) one_sided(z[[2L]], e)[["gap"]], 0)
  cat(sprintf("%-26s %11.3e %11.3e %11.3e %11.3e\n", z[[1L]], g[1],
              g[2], g[3], g[4]))
}
cat("\nA gap that does NOT fall is a corner. A gap that falls needs a\n")
cat("control, because a ONE-SIDED first difference has its own error\n")
cat("of f''(r) * eps, so at a junction with curvature the probe reads\n")
cat("the curvature and cannot separate C1 from C3. Below the cap the\n")
cat("factor is r / (1 - r), whose second derivative is 2 / (1 - r)^3:\n")
cat(sprintf("%-26s %13s %13s\n", "junction", "gap/eps at 1e-5",
            "2 / (1 - r)^3"))
for (z in list(list("r = rlow", RLOW), list("r = rcap", RCAP))) {
  g <- one_sided(z[[2L]], 1e-5)[["gap"]] / 1e-5
  cat(sprintf("%-26s %13.4g %13.4g\n", z[[1L]], g,
              2 / (1 - z[[2L]])^3))
}
cat("\nthe two one-sided derivatives at r = rlow, eps 1e-5:\n")
s <- one_sided(RLOW, 1e-5)
cat(sprintf("  below %.6f  above %.6f  gap %.6g\n", s[["lo"]],
            s[["hi"]], s[["gap"]]))
cat(sprintf("  round 1 read 0.055400 and 1.108040, a gap of %.6f,\n",
            1 / (1 - RLOW)))
cat("  constant over four decades.\n")

cat("\n== does the gate still overshoot the tail? ==\n")
rg <- seq(1e-6, RLOW, length.out = 20001)
f <- vapply(rg, function(r) FF(r)[["f"]], 0)
tail <- rg / (1 - rg)
cat(sprintf("  max f / tail below rlow: %.4f at r = %.5f\n",
            max(f / tail), rg[which.max(f / tail)]))
cat("  a gate may only take the tail away, never add to it, so this\n")
cat("  must be at most 1.\n")
