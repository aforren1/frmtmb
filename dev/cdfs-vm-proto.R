# Prototype: the von Mises tail P(T > u), u in [0, pi], as the angle
# integral exp(k (cos u - 1)) / (2 pi I0e(k)) int_u^hi
# exp(-2 k sin((s + u) / 2) sin((s - u) / 2)) ds by tanh-sinh, the
# range cut where the exponent passes 40.
# Against adaptive quadrature of the density with its peak
# divided out. Output: dev/cdfs-out/vm-proto.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("R/families-cdf.R", envir = E)
ad_pos <- get("ad_pos", E)
ts_rule <- function(h, K) {
  k <- seq(-K, K) * h
  x <- tanh(pi / 2 * sinh(k))
  w <- h * pi / 2 * cosh(k) / cosh(pi / 2 * sinh(k))^2
  list(u = (x + 1) / 2, w = w / 2)
}
vm_tail <- function(u, k, rule) {
  # the exponent 2 k sin((s + u) / 2) sin((s - u) / 2) is at least
  # 2 k sin((s - u) / 2)^2, so it passes 40 within pi sqrt(20 / k) of u
  hi <- pi - ad_pos(pi - u - pi * sqrt(20 / (k + 1e-300)))
  L <- hi - u
  J <- 0
  for (j in seq_along(rule$u)) {
    s <- u + L * rule$u[j]
    J <- J + rule$w[j] * exp(-2 * k * sin((s + u) / 2) * sin((s - u) / 2))
  }
  J <- L * J
  -2 * k * sin(u / 2)^2 - log(2 * pi * besselI(k, 0, expon.scaled = TRUE)) +
    log(1e-300 + J)
}
ref_tail <- function(u, k) {
  if (u >= pi) return(-Inf)
  f <- function(t) exp(k * (cos(t) - cos(u)))
  tot <- stats::integrate(f, u, pi, rel.tol = 1e-13, abs.tol = 0,
                          subdivisions = 5000L)$value
  k * (cos(u) - 1) + log(tot) - log(2 * pi * besselI(k, 0, TRUE))
}
rel <- function(v, r) abs(v - r) / pmax(1, abs(r))
uu <- c(0, 1e-6, 1e-3, 0.01, 0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.1, pi - 1e-3,
        pi - 1e-6)
kk <- c(1e-8, 0.01, 0.3, 1, 3, 9, 10.5, 30, 100, 300, 3000, 3e4)
g <- expand.grid(u = uu, k = kk)
g$ref <- mapply(function(u, k) tryCatch(ref_tail(u, k), error = function(e) NA),
                g$u, g$k)
cat("cases", nrow(g), " without reference", sum(is.na(g$ref)), "\n")
for (cfg in list(list("h0.1", ts_rule(0.1, 30)),
                 list("h0.075", ts_rule(0.075, 40)))) {
  v <- mapply(vm_tail, g$u, g$k, MoreArgs = list(rule = cfg[[2]]))
  e <- rel(v, g$ref)
  i <- which.max(e)
  cat(sprintf("%-7s worst %.2e at u = %g k = %g (ref %.6g, got %.6g); nonfinite %d\n",
              cfg[[1]], max(e, na.rm = TRUE), g$u[i], g$k[i], g$ref[i], v[i],
              sum(!is.finite(v))))
}
