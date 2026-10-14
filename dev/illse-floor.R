# Lane illse: what one double-precision covariance can do on an
# ill-conditioned design, before any frmtmb code. For a raw polynomial
# of degree d on [1, 2] (n = 200, seed 1) and for a covariate far from
# zero, the exact covariance (X'X)^-1 in 256-bit arithmetic (Rmpfr) is
# compared with
#   - lm()'s vcov() (sigma 1), entry by entry;
#   - a' V a for data-range rows a, formed in double from the exact
#     covariance rounded to double, from lm()'s vcov(), and from a
#     factor (a' T)(T' a) with T = R^-1 from lm()'s QR.
# kappa is the condition number of X'X scaled to unit diagonal.
#   Rscript dev/illse-floor.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-mp.R")



rel <- function(a, b) abs(a - b) / abs(b)
one <- function(lab, X, rows) {
  n <- nrow(X)
  XtX <- crossprod(mpfrArray(X, mp_prec, dim = dim(X)))
  Ve <- mp_solve(XtX)
  Vd <- asNumeric(Ve)
  y <- drop(X %*% rep(0, ncol(X))) + stats::rnorm(n)
  fit <- stats::lm.fit(X, y)
  # lm's own vcov at sigma 1: chol2inv of its R
  R <- qr.R(fit$qr)
  piv <- fit$qr$pivot
  if (any(piv != seq_along(piv)) || fit$rank < ncol(X)) {
    cat(sprintf("%-28s lm drops a column (rank %d of %d)\n", lab,
                fit$rank, ncol(X)))
    return(invisible())
  }
  Vlm <- chol2inv(R)
  Dd <- sqrt(diag(asNumeric(XtX)))
  kap <- kappa(asNumeric(XtX) / outer(Dd, Dd), exact = TRUE)
  ent_lm <- mp_entry_err(Vlm, Ve)[["all"]]
  Rinv <- backsolve(R, diag(ncol(X)))
  Rm <- mpfrArray(rows, mp_prec, dim = dim(rows)); qe <- asNumeric(rowSums((Rm %*% Ve) * Rm))
  q_round <- rowSums((rows %*% Vd) * rows)
  q_lm <- rowSums((rows %*% Vlm) * rows)
  q_fac <- rowSums((rows %*% Rinv)^2)
  cat(sprintf(paste0("%-28s kappa %.2e | vcov(lm) entry rel err %.1e |",
                     " a'Va rel err: rounded exact %s; lm %s; factor %s\n"),
              lab, kap, ent_lm,
              format(max(rel(q_round, qe)), digits = 2),
              format(max(rel(q_lm, qe)), digits = 2),
              format(max(rel(q_fac, qe)), digits = 2)))
}
set.seed(1)
x <- stats::runif(200, 1, 2)
xr <- c(1, 1.25, 1.5, 1.75, 2)
for (d in 3:9) {
  X <- outer(x, 0:d, `^`)
  one(sprintf("raw poly degree %d on [1,2]", d), X, outer(xr, 0:d, `^`))
}
for (off in c(1e4, 1e6, 1e7, 1e8, 1e9)) {
  z <- off + stats::runif(200)
  X <- cbind(1, z)
  one(sprintf("x on [%g, %g + 1]", off, off), X,
      cbind(1, off + c(0, 0.25, 0.5, 0.75, 1)))
}
# near-collinear pair x2 = x1 + delta * noise
for (dl in c(1e-4, 1e-6, 1e-8)) {
  x1 <- stats::rnorm(200)
  x2 <- x1 + dl * stats::rnorm(200)
  X <- cbind(1, x1, x2)
  one(sprintf("x2 = x1 + %g noise", dl), X, X[1:5, ])
}
