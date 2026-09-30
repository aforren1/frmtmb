# Can a fixed-depth continued fraction, switched to the complement where
# it does not converge, give log I_x(a, b) to round-off for every
# x < 1/2? Plain doubles against stats::pbeta(log.p = TRUE).
# Seed 20260929; x = kappa / (1 + 2 kappa) with log kappa ~ U(-8, 3),
# log a and log b ~ U(-4, 10).

# log of Lentz's continued fraction for I_x(a, b) (Numerical Recipes
# betacf), a fixed number of steps, the product summed in logs so that
# a branch that does not converge stays finite
cf_log <- function(x, a, b, N) {
  qab <- a + b; qap <- a + 1; qam <- a - 1
  c <- 1
  d <- 1 / (1 - qab * x / qap)
  lh <- log(abs(d))
  for (m in seq_len(N)) {
    m2 <- 2 * m
    aa <- m * (b - m) * x / ((qam + m2) * (a + m2))
    d <- 1 / (1 + aa * d); c <- 1 + aa / c; lh <- lh + log(abs(d * c))
    aa <- -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
    d <- 1 / (1 + aa * d); c <- 1 + aa / c; lh <- lh + log(abs(d * c))
  }
  lh
}
lbeta_i <- function(x, a, b, N) {
  a * log(x) + b * log1p(-x) - lbeta(a, b) - log(a) + cf_log(x, a, b, N)
}
log_ibeta <- function(x, a, b, N) {
  m <- (a + 1) / (a + b + 2)
  t <- m - x
  xd <- 0.5 * (x + m - abs(m - x))
  xc <- 0.5 * (x + m + abs(x - m))
  w <- (t + abs(t)) / (2 * abs(t) + 1e-300)
  ld <- lbeta_i(xd, a, b, N)
  lc <- lbeta_i(1 - xc, b, a, N)
  cap <- log1p(-2^-53)
  lc <- 0.5 * (lc + cap - abs(cap - lc))
  w * ld + (1 - w) * log1p(-exp(lc))
}
set.seed(20260929)
M <- 20000
lk <- runif(M, -8, 3)
x <- exp(lk) / (1 + 2 * exp(lk))
a <- exp(runif(M, -4, 10)); b <- exp(runif(M, -4, 10))
ref <- pbeta(x, a, b, log.p = TRUE)
keep <- is.finite(ref) & ref > -700
cat("points", M, "with -700 < log P:", sum(keep), "\n")
for (N in c(25, 50, 100, 200, 400)) {
  v <- suppressWarnings(log_ibeta(x, a, b, N))
  err <- abs(v - ref) / pmax(1, abs(ref))
  err[!is.finite(err)] <- Inf
  e <- err[keep]
  cat(sprintf("N %4d: rel err max %.2e q99 %.2e q999 %.2e n>1e-12 %5d nonfinite %d\n",
              N, max(e), quantile(e, 0.99), quantile(e, 0.999),
              sum(e > 1e-12), sum(!is.finite(e))))
  if (N %in% c(100, 400)) {
    bad <- keep & err > 1e-12
    if (any(bad)) {
      print(head(data.frame(x = x[bad], a = a[bad], b = b[bad],
                            ref = ref[bad], got = v[bad],
                            err = err[bad])[order(-err[bad]), ], 12))
    }
  }
}
