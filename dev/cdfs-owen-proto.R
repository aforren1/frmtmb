# Prototype quadratures for log G(z, a) = log((1/pi) int_atan(a)^(pi/2)
# exp(-z^2 / (2 cos^2 t)) dt), the skew-normal tail on the side of z.
# Exact references: a = 1 gives F(z; 1) = Phi(z)^2 for z < 0, so
# G(z, 1) = Phi(z)^2 there and G(z, -1) = 1 - Phi(-z)^2 ... ; a = 0
# gives Phi(-|z|); in general adaptive quadrature. Output:
# dev/cdfs-out/owen-proto.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
ns <- asNamespace("frmtmb")
ad_pos <- ns$ad_pos
ad_max <- ns$ad_max
ts_rule <- function(h, K) {
  k <- seq(-K, K) * h
  x <- tanh(pi / 2 * sinh(k))
  w <- h * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2
  list(u = (x + 1) / 2, w = w / 2)
}
gl_unit <- ns$gl_unit
if (is.null(gl_unit)) {
  E <- new.env(parent = ns)
  sys.source("R/families-cdf.R", envir = E)
  gl_unit <- E$gl_unit
}
gl_comp <- function(br, n) {
  g <- gl_unit(n)
  lo <- br[-length(br)]
  wd <- diff(br)
  list(u = as.vector(outer(g$x, wd) + rep(lo, each = n)),
       w = as.vector(outer(g$w, wd)))
}
owen <- function(z, a, rule) {
  ap <- ad_pos(a)
  az <- abs(z)
  z2 <- z * z
  tp <- atan(ap)
  rgt <- pi / 2 - atan(az / sqrt(ap * ap * z2 + 80))
  lft <- ad_max(atan(a), -(pi / 2 - atan(az / sqrt(80))))
  tot <- 0
  for (k in seq_along(rule$u)) {
    t1 <- tp + (rgt - tp) * rule$u[k]
    t2 <- tp - (tp - lft) * rule$u[k]
    g1 <- exp(-0.5 * z2 * (tan(t1)^2 - ap * ap))
    g2 <- exp(-0.5 * z2 * (tan(t2)^2 - ap * ap))
    tot <- tot + rule$w[k] * ((rgt - tp) * g1 + (tp - lft) * g2)
  }
  -log(pi) - 0.5 * z2 * (1 + ap * ap) + log(tot)
}
ref_G <- function(z, a) {
  # adaptive Gauss-Kronrod on the angle, the peak divided out, split at
  # the peak; independent of the fixed rules it checks
  ap <- max(a, 0)
  f <- function(t) exp(-z^2 / 2 * (tan(t)^2 - ap^2))
  tp <- atan(ap)
  i1 <- stats::integrate(f, tp, pi / 2, rel.tol = 1e-13, abs.tol = 0,
                         subdivisions = 2000L)$value
  i2 <- if (a < 0) stats::integrate(f, atan(a), 0, rel.tol = 1e-13,
                                    abs.tol = 0, subdivisions = 2000L)$value
  else 0
  log(i1 + i2) - log(pi) - z^2 * (1 + ap^2) / 2
}
rules <- list(
  geo6x10 = gl_comp(c(0, 1/32, 1/16, 1/8, 1/4, 1/2, 1), 10),
  geo2side = gl_comp(c(0, 1/32, 1/16, 1/8, 1/4, 1/2, 3/4, 7/8, 15/16,
                       31/32, 1), 8),
  ts_h0.1 = ts_rule(0.1, 30),
  ts_h0.075 = ts_rule(0.075, 40),
  ts_h0.05 = ts_rule(0.05, 62)
)
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
cases <- expand.grid(z = c(-38, -20, -8, -3, -1, -0.3, -0.03, -1e-4, 0, 1e-4,
                           0.03, 0.3, 1, 3, 8, 20, 38),
                     a = c(-30, -3, -1, -0.2, 0, 0.2, 1, 3, 30))
# exact where the parameter makes it closed form
exact <- function(z, a) {
  if (a == 0) return(pnorm(-abs(z), log.p = TRUE))
  # G(z, 1) is F(z; 1) = Phi(z)^2 for z < 0 and S(z; -1) = Phi(-z)^2 for
  # z > 0: Phi(-|z|)^2 either way
  if (a == 1) return(2 * pnorm(-abs(z), log.p = TRUE))
  # G(z, -1) = F(z; -1) for z < 0 = 1 - Phi(-z)^2 = Phi(z) (1 + Phi(-z))
  if (a == -1) return(pnorm(-abs(z), log.p = TRUE) + log1p(pnorm(abs(z))))
  NA
}
cases$ref <- mapply(function(z, a) {
  e <- exact(z, a)
  if (!is.na(e)) e else ref_G(z, a)
}, cases$z, cases$a)
for (nm in names(rules)) {
  v <- mapply(owen, cases$z, cases$a, MoreArgs = list(rule = rules[[nm]]))
  e <- rel(v, cases$ref)
  i <- which.max(e)
  cat(sprintf("%-10s nodes %3d  worst %.2e at z = %g a = %g; exact-a cases worst %.2e\n",
              nm, 2 * length(rules[[nm]]$u), max(e), cases$z[i], cases$a[i],
              max(e[cases$a %in% c(-1, 0, 1)])))
}
