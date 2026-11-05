# Reads, at 256 bits, the exact point the double-precision closed form calls
# its maximum on y = exp(0.8 x) + 0.3 e, seeds 5701 and 5702, and a line search
# by exact arithmetic around frm()'s stop: is the 0.68 nats real?
#   Rscript dev/ncgp-exp8-mpfr2.R <seed>
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(Rmpfr))
seed <- as.integer(commandArgs(TRUE)[1])
set.seed(seed)
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
D2 <- outer(x, x, "-")^2
marg <- function(p) {
  s2 <- exp(2 * p[2]); K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + 1e-12)
  diag(K) <- diag(K) + exp(2 * p[4])
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
start <- c(mean(y), log(5e4), log(4.1), log(0.33))
sc <- c(sd(y), 1, 1, 1)
fn <- function(q) marg(q * sc)
o <- nlminb(start / sc, fn, control = list(iter.max = 800, eval.max = 1500, rel.tol = 1e-12))
o2 <- optim(o$par, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-14))
best <- if (o2$value < o$objective) o2$par else o$par
pb <- best * sc
cat(sprintf("NCGPX2 seed %d double best ll %.6f at %s\n", seed, -min(o2$value, o$objective),
            paste(sprintf("%.12g", pb), collapse = " ")))
cat(sprintf("NCGPX2 exact at that point %.6f\n", -exact(pb)))
for (h in c(1e-8, 1e-6, 1e-4)) {
  pp <- pb; pp[3] <- pp[3] + h
  cat(sprintf("NCGPX2 lscale + %.0e: double %.6f exact %.6f\n", h, -marg(pp), -exact(pp)))
}
