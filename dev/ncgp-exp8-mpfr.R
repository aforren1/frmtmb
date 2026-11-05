# Is the closed-form maximum on y = exp(0.8 x) + 0.3 e (seed 5701) real, or the
# rounding of a kernel conditioned 2e12? Reads the marginal at 256 bits at the
# frm() stop, at the closed form's best, and on the line between. Rmpfr.
#   Rscript dev/ncgp-exp8-mpfr.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(Rmpfr))
set.seed(5701)
n <- 100
x <- sort(runif(n, 0, 10))
e0 <- rnorm(n)
y <- exp(0.8 * x) + 0.3 * e0
prec <- 256
mx <- mpfr(x, prec); my <- mpfr(y, prec)
exact <- function(p) {
  s2 <- exp(2 * mpfr(p[2], prec)); r2 <- exp(2 * mpfr(p[3], prec))
  sg2 <- exp(2 * mpfr(p[4], prec))
  V <- array(mpfr(0, prec), c(n, n))
  for (i in seq_len(n)) for (j in seq_len(i)) {
    v <- s2 * exp(-(mx[i] - mx[j])^2 / (2 * r2))
    if (i == j) v <- v * (1 + mpfr(1e-12, prec)) + sg2
    V[i, j] <- v
  }
  # LDL, exact log det and quadratic form
  L <- array(mpfr(0, prec), c(n, n))
  r <- my - mpfr(p[1], prec)
  z <- mpfr(rep(0, n), prec)
  ld <- mpfr(0, prec)
  for (j in seq_len(n)) {
    s <- V[j, j]
    if (j > 1) for (k in seq_len(j - 1L)) s <- s - L[j, k]^2
    L[j, j] <- sqrt(s)
    ld <- ld + log(L[j, j])
    if (j < n) for (i in (j + 1L):n) {
      t <- V[i, j]
      if (j > 1) for (k in seq_len(j - 1L)) t <- t - L[i, k] * L[j, k]
      L[i, j] <- t / L[j, j]
    }
    t <- r[j]
    if (j > 1) for (k in seq_len(j - 1L)) t <- t - L[j, k] * z[k]
    z[j] <- t / L[j, j]
  }
  as.numeric(0.5 * sum(z^2) + ld + 0.5 * n * log(2 * pi))
}
marg <- function(p) {
  D2 <- outer(x, x, "-")^2
  s2 <- exp(2 * p[2]); K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + 1e-12)
  diag(K) <- diag(K) + exp(2 * p[4])
  R <- chol(K)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
pf <- c(364.760695, 10.885579, 1.402627, -1.079025)
pc <- c(364.760695, log(53373), log(4.1354), log(0.3408))
for (a in c(0, 0.25, 0.5, 0.75, 1)) {
  p <- pf + a * (pc - pf)
  cat(sprintf("NCGPX a=%.2f double %.6f exact %.6f\n", a, -marg(p), -exact(p)))
}
