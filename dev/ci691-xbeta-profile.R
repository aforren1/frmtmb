# Lane ci691: the profile log-likelihood of the test's xbeta data
# (test-xbeta-zibb-hurdle-cum.R, seed 4200, phi 1e4 and 2e4) along
# log kappa, from the stats-function reference alone, to see whether the
# maximum is interior or at the normal limit kappa -> Inf.
#   Rscript dev/ci691-xbeta-profile.R
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi
  b <- (1 - mu) * phi
  d <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / d, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}
ll <- function(p, y) sum(ref_xbeta(y, plogis(p[1]), exp(p[2]), exp(p[3])))
for (phi in c(1e4, 2e4)) {
  set.seed(4200)
  z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
  y <- pmin(pmax(1.1 * z - 0.05, 0), 1)
  cat("\n==== phi", phi, " zeros", sum(y == 0), " ones", sum(y == 1), "\n")
  prev <- c(qlogis(0.047), log(phi))
  for (lk in seq(-4, 3, by = 0.25)) {
    # carry the previous optimum along the ridge as the start
    f <- function(q) -ll(c(q, lk), y)
    o <- nlminb(prev, f, control = list(eval.max = 5000, iter.max = 5000,
                                       rel.tol = 1e-14))
    prev <- o$par
    cat(sprintf("log kappa %6.2f  logit mu %9.5f  log phi %9.5f  ll %.8f\n",
                lk, o$par[1], o$par[2], -o$objective))
  }
}
