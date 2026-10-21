# Prototype 2 for log G(z, a) = log(Phi(-|z|) - 2 T(z, a)), the
# skew-normal tail on the side of z. Two regimes:
#  - |z| <= 1.2: Owen's T by a 20-point Gauss-Legendre rule in x on
#    [0, a] for |a| <= 1.1, and through T(h, a) = Phi(h)/2 + Phi(ah)/2 -
#    Phi(h) Phi(ah) - T(ah, 1/a) above; G is O(1) there, so the
#    subtraction costs nothing relative;
#  - |z| >= 0.8: the angle integral with the peak divided out, by a
#    composite rule geometric towards the peak, or tanh-sinh.
# Exact references at a = 0, +-1; elsewhere stats::integrate() on the
# angle at rel.tol 1e-13. Output: dev/cdfs-out/owen-proto2.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("R/families-cdf.R", envir = E)
ad_pos <- get("ad_pos", E)
ad_max <- get("ad_max", E)
gl_unit <- get("gl_unit", E)
gl_comp <- function(br, n) {
  g <- gl_unit(n)
  lo <- br[-length(br)]
  wd <- diff(br)
  list(u = as.vector(outer(g$x, wd) + rep(lo, each = n)),
       w = as.vector(outer(g$w, wd)))
}
ts_rule <- function(h, K) {
  k <- seq(-K, K) * h
  x <- tanh(pi / 2 * sinh(k))
  w <- h * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2
  list(u = (x + 1) / 2, w = w / 2)
}
owen_theta <- function(z, a, rule) {
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
g20 <- gl_unit(20L)
owen_T_gl <- function(h, a) {
  # (1 / 2 pi) int_0^a exp(-h^2 (1 + x^2) / 2) / (1 + x^2) dx
  tot <- 0
  for (k in seq_along(g20$x)) {
    x <- a * g20$x[k]
    tot <- tot + g20$w[k] * exp(-0.5 * h * h * (1 + x * x)) / (1 + x * x)
  }
  a * tot / (2 * pi)
}
owen_T <- function(h, a) {
  if (abs(a) <= 1.1) return(owen_T_gl(h, a))
  b <- abs(a)
  ph <- pnorm(h)
  pah <- pnorm(b * h)
  sign(a) * (0.5 * ph + 0.5 * pah - ph * pah - owen_T_gl(b * h, 1 / b))
}
owen_small <- function(z, a) log(pnorm(-abs(z)) - 2 * owen_T(z, a))
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
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
av <- c(-1e3, -30, -3, -1.1, -1, -0.95, -0.2, 0, 0.2, 0.95, 1, 1.1, 3, 30,
        1e3)
ref_of <- function(z, a) {
  e <- exact(z, a)
  if (!is.na(e)) return(e)
  tryCatch(ref_G(z, a), error = function(err) NA_real_)
}
small <- expand.grid(z = c(-1.2, -0.8, -0.3, -0.03, -1e-4, 0, 1e-4, 0.03,
                           0.3, 0.8, 1.2), a = av)
small$ref <- mapply(ref_of, small$z, small$a)
v <- mapply(owen_small, small$z, small$a)
e <- rel(v, small$ref)
cat(sprintf("small |z|: %d cases, %d without a reference; worst %.2e at z = %g a = %g\n",
            nrow(small), sum(is.na(small$ref)), max(e, na.rm = TRUE),
            small$z[which.max(e)], small$a[which.max(e)]))
big <- expand.grid(z = c(-38, -25, -12, -6, -3, -1.5, -1, -0.8, 0.8, 1, 1.5,
                         3, 6, 12, 25, 38), a = av)
big$ref <- mapply(ref_of, big$z, big$a)
rules <- list(
  geo6x10 = gl_comp(c(0, 1/32, 1/16, 1/8, 1/4, 1/2, 1), 10),
  geo6x8 = gl_comp(c(0, 1/32, 1/16, 1/8, 1/4, 1/2, 1), 8),
  geo8x8 = gl_comp(c(0, 1/128, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1), 8),
  ts_h0.1 = ts_rule(0.1, 30),
  ts_h0.075 = ts_rule(0.075, 40)
)
for (nm in names(rules)) {
  v <- mapply(owen_theta, big$z, big$a, MoreArgs = list(rule = rules[[nm]]))
  e <- rel(v, big$ref)
  ex <- big$a %in% c(-1, 0, 1)
  cat(sprintf("|z| >= 0.8 %-9s nodes %3d: %d cases, worst %.2e at z = %g a = %g; exact-a worst %.2e; nonfinite %d\n",
              nm, 2 * length(rules[[nm]]$u), sum(!is.na(big$ref)),
              max(e, na.rm = TRUE), big$z[which.max(e)], big$a[which.max(e)],
              max(e[ex]), sum(!is.finite(v))))
}
