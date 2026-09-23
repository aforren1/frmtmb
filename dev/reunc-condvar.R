# Lane wt-reunc: what the last round's "add Var(b | y)" arms measured.
#
# dev/shapes-coverage.R (condvar, 0.9609) and dev/shapes-rev-cov.R
# (joint_plus_condvar, 0.9477) both widened the modes-conditional
# predict() interval by v_g = sigma^2 tau^2 / (sigma^2 + n_g tau^2), the
# conditional variance of u_g GIVEN beta. The interval they widened
# already carried the fixed-effect variance x' Var(beta_hat) x. The
# exact prediction error variance is not the sum of the two: the BLUP
# b_g = lambda_g (ybar_g - xbar_g' beta_hat) moves AGAINST the
# intercept, so the error is (x - lambda_g xbar_g)' (beta_hat - beta)
# plus the BLUP error at known beta. Adding v_g to x' V x counts the
# intercept's share twice.
#
# This computes, at the TRUE variance parameters and exactly (Henderson's
# mixed model equations), the coverage each construction would have if
# the variance parameters were known: the "sum" interval
# sigma^2 + x'Vx + v_g against the exact sigma^2 + PEV. Averaged over
# the designs' own data and new-point distributions, 4000 replicates.
#
#   Rscript dev/reunc-condvar.R
cover_of <- function(design, nrep = 4000L, seed = 1L) {
  set.seed(seed)
  p <- switch(design,
    shapes = list(n = 60, ng = 12, sigma = 1, tau = 0.7),
    review = list(n = 60, ng = 12, sigma = 0.8, tau = 0.7),
    few4 = list(n = 40, ng = 4, sigma = 1, tau = 0.7))
  z <- stats::qnorm(0.975)
  cov_sum <- numeric(0)
  ratio <- numeric(0)
  for (r in seq_len(nrep)) {
    x <- stats::rnorm(p$n)
    g <- factor(rep(seq_len(p$ng), length.out = p$n))
    X <- cbind(1, x)
    Z <- stats::model.matrix(~ 0 + g)
    s2 <- p$sigma^2
    t2 <- p$tau^2
    C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
               cbind(crossprod(Z, X), crossprod(Z) + diag(s2 / t2, p$ng))) /
      s2
    Ci <- solve(C)
    # Var(beta_hat) of GLS at known variances is the beta block of C^-1
    Vb <- Ci[1:2, 1:2]
    gsel <- sample(seq_len(p$ng), 10, TRUE)
    xn <- stats::rnorm(10)
    An <- cbind(1, xn, diag(p$ng)[gsel, , drop = FALSE])
    pev <- rowSums((An %*% Ci) * An)
    ng_obs <- as.numeric(table(g))[gsel]
    v <- s2 * t2 / (s2 + ng_obs * t2)
    xv <- rowSums((cbind(1, xn) %*% Vb) * cbind(1, xn))
    w_sum <- sqrt(s2 + xv + v)
    w_exact <- sqrt(s2 + pev)
    # the error is exactly N(0, s2 + pev), so an interval of half-width
    # z * w covers with probability 2 Phi(z w / w_exact) - 1
    cov_sum <- c(cov_sum, 2 * stats::pnorm(z * w_sum / w_exact) - 1)
    ratio <- c(ratio, (xv + v) / pev)
  }
  cat(sprintf(paste0("%-7s sigma %.1f tau %.1f groups %2d: known-parameter ",
                     "coverage of sigma^2 + x'Vx + v_g = %.4f (exact PEV ",
                     "interval 0.9500); (x'Vx + v_g) / PEV mean %.3f\n"),
              design, p$sigma, p$tau, p$ng, mean(cov_sum), mean(ratio)))
}
cat("script dev/reunc-condvar.R, seed 1, 4000 replicates of 10 points\n")
for (d in c("shapes", "review", "few4")) cover_of(d)
