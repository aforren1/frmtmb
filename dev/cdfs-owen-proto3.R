# Prototype 3: the combined log G(z, a). The regime is chosen by the
# tail depth t = |z| sqrt(1 + max(a, 0)^2): Owen's T by Gauss-Legendre
# where t is small (G is then at least 0.1, so absolute accuracy is
# relative accuracy), tanh-sinh on the angle where it is not. Grid
# includes small |z| with large a, where t is large. Output:
# dev/cdfs-out/owen-proto3.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("R/families-cdf.R", envir = E)
ad_pos <- get("ad_pos", E)
ad_max <- get("ad_max", E)
ad_min <- get("ad_min", E)
ad_smooth01 <- get("ad_smooth01", E)
gl_unit <- get("gl_unit", E)
ts_rule <- function(h, K) {
  k <- seq(-K, K) * h
  x <- tanh(pi / 2 * sinh(k))
  w <- h * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2
  list(u = (x + 1) / 2, w = w / 2)
}
g20 <- gl_unit(20L)
owen_T_gl <- function(h, a) {
  tot <- 0
  for (k in seq_along(g20$x)) {
    x <- a * g20$x[k]
    tot <- tot + g20$w[k] * exp(-0.5 * h * h * (1 + x * x)) / (1 + x * x)
  }
  a * tot / (2 * pi)
}
owen_T_ad <- function(h, a) {
  aa <- abs(a)
  wI <- ad_smooth01((aa - 0.9) / 0.2)
  ad <- ad_max(ad_min(a, 1.1), -1.1)
  b <- ad_max(aa, 0.9)
  ph <- RTMB::pnorm(h)
  pbh <- RTMB::pnorm(b * h)
  ident <- sign(a) * (0.5 * ph + 0.5 * pbh - ph * pbh - owen_T_gl(b * h, 1 / b))
  (1 - wI) * owen_T_gl(h, ad) + wI * ident
}
owen_ts <- function(z, a, rule) {
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
owen_comb <- function(z, a, rule, t0, t1) {
  ap <- ad_pos(a)
  sq <- sqrt(1 + ap * ap)
  t <- abs(z) * sq
  ws <- ad_smooth01((t1 - t) / (t1 - t0))
  zs <- sign(z) * ad_min(abs(z), t1 / sq)
  small <- log(RTMB::pnorm(-abs(zs)) - 2 * owen_T_ad(zs, a))
  ws * small + (1 - ws) * owen_ts(z, a, rule)
}
ref_G <- function(z, a) {
  ap <- max(a, 0)
  f <- function(t) exp(-z^2 / 2 * (tan(t)^2 - ap^2))
  tp <- atan(ap)
  i1 <- stats::integrate(f, tp, pi / 2, rel.tol = 1e-13, abs.tol = 0,
                         subdivisions = 5000L)$value
  i2 <- if (a < 0) stats::integrate(f, atan(a), 0, rel.tol = 1e-13,
                                    abs.tol = 0, subdivisions = 5000L)$value
  else 0
  log(i1 + i2) - log(pi) - z^2 * (1 + ap^2) / 2
}
exact <- function(z, a) {
  if (a == 0) return(pnorm(-abs(z), log.p = TRUE))
  if (a == 1) return(2 * pnorm(-abs(z), log.p = TRUE))
  if (a == -1) return(pnorm(-abs(z), log.p = TRUE) + log1p(pnorm(abs(z))))
  NA
}
ref_of <- function(z, a) {
  e <- exact(z, a)
  if (!is.na(e)) return(e)
  tryCatch(ref_G(z, a), error = function(err) NA_real_)
}
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
av <- c(-1e3, -30, -3, -1.1, -1, -0.95, -0.2, 0, 0.2, 0.95, 1, 1.1, 3, 30,
        1e3)
zv <- c(-38, -12, -3, -1.5, -1, -0.5, -0.1, -0.01, -1e-4, 0, 1e-4, 0.01,
        0.1, 0.5, 1, 1.5, 3, 12, 38)
g <- expand.grid(z = zv, a = av)
g$ref <- mapply(ref_of, g$z, g$a)
g$t <- abs(g$z) * sqrt(1 + pmax(g$a, 0)^2)
cat("cases", nrow(g), " without reference", sum(is.na(g$ref)), "\n")
for (cfg in list(list("ts_h0.1", ts_rule(0.1, 30), 0.8, 1.2),
                 list("ts_h0.1", ts_rule(0.1, 30), 1.2, 1.6),
                 list("ts_h0.075", ts_rule(0.075, 40), 0.8, 1.2),
                 list("ts_h0.075", ts_rule(0.075, 40), 1.2, 1.6))) {
  v <- mapply(owen_comb, g$z, g$a,
              MoreArgs = list(rule = cfg[[2]], t0 = cfg[[3]], t1 = cfg[[4]]))
  e <- rel(v, g$ref)
  i <- which.max(e)
  cat(sprintf("%-9s band [%g, %g]: worst %.2e at z = %g a = %g (t = %.3g); exact-a worst %.2e; nonfinite %d\n",
              cfg[[1]], cfg[[3]], cfg[[4]], max(e, na.rm = TRUE), g$z[i],
              g$a[i], g$t[i], max(e[g$a %in% c(-1, 0, 1)]), sum(!is.finite(v))))
}
