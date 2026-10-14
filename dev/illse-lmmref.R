# Lane illse: the reference lmm_exact_cov() (dev/illse-common.R, its
# sufficient-statistic form) against the dense marginal likelihood in
# base R double with optimHess() on a small well-conditioned design
# (n 40, 4 groups, seed 1): the two must agree to the finite-difference
# Hessian's own accuracy (about 1e-7 relative).
#   Rscript dev/illse-lmmref.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
set.seed(1)
n <- 40
g <- factor(rep(1:4, each = 10))
X <- cbind(1, rnorm(n))
y <- drop(X %*% c(1, 2)) + rnorm(4)[g] + rnorm(n)
par <- c(1.1, 1.9, log(0.9), log(0.8))
Z <- model.matrix(~ g - 1)
nll <- function(th) {
  V <- exp(2 * th[3]) * diag(n) + exp(2 * th[4]) * tcrossprod(Z)
  r <- y - X %*% th[1:2]
  0.5 * (determinant(V)$modulus + sum(r * solve(V, r)))
}
Vd <- solve(optimHess(par, nll))
Vm <- asNumeric(lmm_exact_cov(X, y, g, par))
cat(sprintf("largest relative difference over the entries: %.2e\n",
            max(abs(Vd - Vm) / abs(Vm))))
# and against the same marginal likelihood written per row in 256 bits
# (the form the sufficient statistics replace), which must agree to far
# below double precision
Xm <- mpfrArray(X, mp_prec, dim = dim(X))
ym <- mp(y)
gi <- split(seq_along(y), g)
nll_row <- function(th) {
  s2 <- exp(2 * th[3])
  t2 <- exp(2 * th[4])
  r <- ym - Xm %*% th[1:2]
  out <- mp(0)
  for (ix in gi) {
    den <- s2 + length(ix) * t2
    out <- out + (length(ix) - 1) * log(s2) + log(den) +
      (sum(r[ix]^2) - t2 / den * sum(r[ix])^2) / s2
  }
  out / 2
}
Vr <- mp_solve(mp_hess(nll_row, par))
cat(sprintf("against the per-row form in 256 bits: %.2e\n",
            max(asNumeric(abs(Vr - lmm_exact_cov(X, y, g, par)) / abs(Vr)))))
