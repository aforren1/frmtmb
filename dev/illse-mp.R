# Lane illse: extended-precision helpers (Rmpfr). Rmpfr has no solve()
# for an mpfrMatrix on this build, so a Gauss-Jordan inverse with
# partial pivoting, adequate for the p <= 30 matrices used here.
suppressMessages(library(Rmpfr))
mp_prec <- 256
mp <- function(x) mpfr(x, mp_prec)
mp_solve <- function(A) {
  if (!inherits(A, "mpfrMatrix")) A <- mpfrArray(A, mp_prec, dim = dim(A))
  p <- nrow(A)
  B <- mpfrArray(diag(p), mp_prec, dim = c(p, p))
  for (k in seq_len(p)) {
    piv <- k - 1L + which.max(asNumeric(abs(A[k:p, k])))
    if (piv != k) {
      tA <- A[k, ]; A[k, ] <- A[piv, ]; A[piv, ] <- tA
      tB <- B[k, ]; B[k, ] <- B[piv, ]; B[piv, ] <- tB
    }
    s <- A[k, k]
    A[k, ] <- A[k, ] / s
    B[k, ] <- B[k, ] / s
    for (i in seq_len(p)[-k]) {
      f <- A[i, k]
      if (asNumeric(f) != 0) {
        A[i, ] <- A[i, ] - f * A[k, ]
        B[i, ] <- B[i, ] - f * B[k, ]
      }
    }
  }
  B
}
# Hessian of an mpfr-capable function by central second differences in
# extended precision: with 256 bits and a step of 1e-20 the truncation
# (h^2) and the rounding (2^-256 / h^2) are both far below double's.
mp_hess <- function(f, x, h = 1e-20) {
  p <- length(x)
  x <- mp(x)
  h <- mp(h)
  H <- mpfrArray(0, mp_prec, dim = c(p, p))
  f0 <- f(x)
  e <- function(i) { v <- mp(rep(0, p)); v[i] <- h; v }
  fp <- lapply(seq_len(p), function(i) f(x + e(i)))
  fm <- lapply(seq_len(p), function(i) f(x - e(i)))
  for (i in seq_len(p)) {
    H[i, i] <- (fp[[i]] - 2 * f0 + fm[[i]]) / h^2
    if (i < p) for (j in (i + 1):p) {
      v <- (f(x + e(i) + e(j)) - f(x + e(i) - e(j)) -
              f(x - e(i) + e(j)) + f(x - e(i) - e(j))) / (4 * h^2)
      H[i, j] <- v
      H[j, i] <- v
    }
  }
  H
}
# largest relative error over the entries of a covariance against an
# mpfr reference, and over the diagonal alone
mp_entry_err <- function(V, Ve) {
  E <- asNumeric(abs(mpfrArray(V, mp_prec, dim = dim(V)) - Ve) / abs(Ve))
  c(all = max(E), diag = max(diag(E)))
}
# errors of a covariance V against an mpfr reference Ve: `cor` the
# largest entry error over sqrt(Ve_ii Ve_jj) (every entry, on the scale
# of the correlation), `rel` the largest relative error over the entries
# whose reference correlation is at least 1e-3 in size, `diag` over the
# variances
mp_cov_err <- function(V, Ve) {
  Vm <- mpfrArray(V, mp_prec, dim = dim(V))
  d <- sqrt(diag(asNumeric(Ve)))
  Ab <- asNumeric(abs(Vm - Ve))
  ve <- asNumeric(Ve)
  cr <- abs(ve) / outer(d, d)
  big <- cr >= 1e-3
  c(cor = max(Ab / outer(d, d)), rel = max((Ab / abs(ve))[big]),
    diag = max(diag(Ab) / diag(abs(ve))))
}
