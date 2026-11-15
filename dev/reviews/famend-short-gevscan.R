# Short check of lane famend, F1 attack: map the outside GEV rows whose
# gradient is still non-finite on the lane build. Where cl - z lies in
# [8, 16) its rounding puts z + out at -1 + 2^-49, 7 ulps ABOVE the
# clamp -1 + 1e-15, so the ratio is not held by gev_ratio_top. The
# test's setting: mu 0, sigma 1, xi through log1p.
#   Rscript dev/reviews/famend-short-gevscan.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
fam <- gen_extreme_value()
grad <- function(y, xi, mu = 0, ls = 0) {
  p <- c(mu, ls, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                     .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  tp$jacobian(p)
}
cl <- -1 + 1e-15
d <- c(8, 16)
for (dd in d) {
  z <- cl - dd
  out <- 0.5 * ((cl - z) + abs(cl - z))
  cat(sprintf("cl - z = %g: z + out - cl = %.4g\n", dd, (z + out) - cl))
}
xs <- seq(-19.6, -20.8, by = -0.01)
tot <- 0
for (xi in seq(0.40, 0.90, by = 0.02)) {
  nb <- 0
  xr <- c()
  for (x in xs) {
    b <- any(!is.finite(grad(x, xi)))
    nb <- nb + b
    if (b) xr <- c(xr, x)
  }
  tot <- tot + nb
  cat(sprintf("xi %.2f: non-finite gradient %3d of %d rows%s\n", xi, nb,
              length(xs), if (nb) sprintf(", x from %.2f to %.2f",
                                          max(xr), min(xr)) else ""))
}
cat("total non-finite rows:", tot, "\n")
# the same band with another mu and sigma (rows in scale units)
nb <- 0
for (x in seq(-19.8, -20.6, by = -0.02)) {
  nb <- nb + any(!is.finite(grad(3 + exp(0.7) * x, 0.6, mu = 3, ls = 0.7)))
}
cat("xi 0.6, mu 3, sigma exp(0.7): non-finite rows", nb, "of 41\n")
