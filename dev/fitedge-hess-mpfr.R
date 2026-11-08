# Item 3: which outer Hessian is right where the two routes disagree? The exact
# Hessian of the gaussian marginal, written from the model statement with no
# frmtmb code, y ~ N(b0, sd^2 (K + 1e-12 I) + sigma^2 I), in 200-bit
# arithmetic (Rmpfr), by central second differences of the value at a step of
# 1e-20. Design: test-gp-start.R's "y's and x's units", seed 1, n 80,
# 2 sin(3x) + x, noise 1e-3, at scale y 1 and y 1e-3. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
suppressPackageStartupMessages(library(Rmpfr))
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
prec <- 200L
nll_mp <- function(p, x, y) {
  n <- length(y)
  b0 <- p[1]; sig2 <- exp(2 * p[2]); sd2 <- exp(2 * p[3]); rho2 <- exp(2 * p[4])
  D2 <- outer(x, x, "-")^2
  A <- mpfrArray(0, prec, dim = c(n, n))
  K <- exp(-mpfr(D2, prec) / (2 * rho2))
  A <- sd2 * K
  diag(A) <- diag(A) + sd2 * mpfr(1e-12, prec) + sig2
  L <- mpfrArray(0, prec, dim = c(n, n))
  for (j in seq_len(n)) {
    s <- A[j, j]
    if (j > 1) s <- s - sum(L[j, 1:(j - 1)]^2)
    L[j, j] <- sqrt(s)
    if (j < n) {
      v <- A[(j + 1):n, j]
      if (j > 1) {
        v <- v - as.vector(L[(j + 1):n, 1:(j - 1), drop = FALSE] %*% L[j, 1:(j - 1)])
      }
      L[(j + 1):n, j] <- v / L[j, j]
    }
  }
  r <- mpfr(y, prec) - b0
  z <- mpfrArray(0, prec, dim = n)
  for (i in seq_len(n)) {
    s <- r[i]
    if (i > 1) s <- s - sum(L[i, 1:(i - 1)] * z[1:(i - 1)])
    z[i] <- s / L[i, i]
  }
  sum(log(diag(L))) + sum(z^2) / 2 + n * log(2 * Const("pi", prec)) / 2
}
hess_mp <- function(p0, x, y) {
  k <- length(p0)
  h <- mpfr(1e-20, prec)
  p <- mpfr(p0, prec)
  f0 <- nll_mp(p, x, y)
  H <- matrix(0, k, k)
  e <- function(i) { v <- mpfr(rep(0, k), prec); v[i] <- h; v }
  fp <- lapply(seq_len(k), function(i) nll_mp(p + e(i), x, y))
  fm <- lapply(seq_len(k), function(i) nll_mp(p - e(i), x, y))
  for (i in seq_len(k)) {
    H[i, i] <- asNumeric((fp[[i]] - 2 * f0 + fm[[i]]) / h^2)
    for (j in seq_len(i - 1)) {
      v <- (nll_mp(p + e(i) + e(j), x, y) - nll_mp(p + e(i) - e(j), x, y) -
              nll_mp(p - e(i) + e(j), x, y) + nll_mp(p - e(i) - e(j), x, y)) /
        (4 * h^2)
      H[i, j] <- H[j, i] <- asNumeric(v)
    }
  }
  list(H = H, f0 = asNumeric(f0))
}
set.seed(1)
x <- round(stats::runif(80L, 0, 10), 3)
y <- 2 * sin(3 * x) + x + 1e-3 * stats::rnorm(80L)
for (ys in c(1, 1e-3)) {
  f <- suppressWarnings(frm(y ~ gp(x), data = data.frame(x = x, y = y * ys)))
  at <- f$opt$par
  env <- f$obj$env
  hz <- {
    unlockBinding("gp_marginal_outer_gr", ns)
    assign("gp_marginal_outer_gr", function(fit, at) NULL, envir = ns)
    lockBinding("gp_marginal_outer_gr", ns)
    rm(list = intersect(".frmtmb_outer_h", ls(env, all.names = TRUE)), envir = env)
    ns$fit_outer_hessian(f)$H
  }
  unlockBinding("gp_marginal_outer_gr", ns)
  assign("gp_marginal_outer_gr", orig, envir = ns)
  lockBinding("gp_marginal_outer_gr", ns)
  # the marginal route whatever the value gate says
  gm <- local({
    tl <- ns$gp_marginal_value_tol
    unlockBinding("gp_marginal_value_tol", ns)
    assign("gp_marginal_value_tol", Inf, envir = ns)
    on.exit({
      assign("gp_marginal_value_tol", tl, envir = ns)
      lockBinding("gp_marginal_value_tol", ns)
    })
    rm(list = intersect(".frmtmb_outer_h", ls(env, all.names = TRUE)), envir = env)
    ns$fit_outer_hessian(f)$H
  })
  t0 <- proc.time()[[3]]
  ex <- hess_mp(as.numeric(at), x, y * ys)
  se <- function(H) sqrt(diag(solve(H)))
  cat(sprintf("MP y scale %g: mpfr %.1f s, value mpfr %.10f fit %.10f\n", ys,
              proc.time()[[3]] - t0, ex$f0, f$opt$objective))
  cat("MP  names", names(at), "\n")
  cat("MP  SE exact   ", signif(se(ex$H), 7), "\n")
  cat("MP  SE z tape  ", signif(se(hz), 7), "\n")
  cat("MP  SE marginal", signif(se(gm), 7), "\n")
  cat(sprintf("MP  max rel SE error: z tape %.3e, marginal %.3e\n",
              max(abs(se(hz) / se(ex$H) - 1)), max(abs(se(gm) / se(ex$H) - 1))))
  cat(sprintf("MP  max |H - H_exact| / max|H_exact|: z tape %.3e, marginal %.3e\n",
              max(abs(hz - ex$H)) / max(abs(ex$H)), max(abs(gm - ex$H)) / max(abs(ex$H))))
}
cat("DONE\n")
