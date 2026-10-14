# Reviewer round 2: the 67 G1 rows (brms's exported d/p/q/r functions)
# are now "no counterpart (decision)", out of the rate's denominator.
# Can the d and r rows for families frmtmb HAS be answered by the same
# thin read the lane uses for G2's dcumulative() (a family object's
# lpdf / the shared simulator)? Each assertion is brms 2.23.0's
# tests.distributions.R text, with brms's function replaced by a shim
# over frmtmb::row_lpdf() or frmtmb::sim_draw(). Seed 20261013.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20261013)
d_ <- function(fun, x, dp, av = list(), log = FALSE) {
  fam <- frmtmb::brmsfamily(fun)
  n <- length(x)
  dp <- lapply(dp, rep_len, length.out = n)
  av <- lapply(av, rep_len, length.out = n)
  ll <- as.numeric(frmtmb::row_lpdf(fam, x, x, dp, av, NULL))
  if (log) ll else exp(ll)
}
r_ <- function(fun, n, dp, av = list()) {
  fam <- frmtmb::brmsfamily(fun)
  dp <- lapply(dp, rep_len, length.out = n)
  av <- lapply(av, rep_len, length.out = n)
  as.vector(frmtmb::sim_draw(list(family = fam, dpars = dp, aterms = av,
                                  n = n, autocor = NULL, extra = NULL,
                                  max_iter = NULL)))
}
res <- list()
chk <- function(id, expr) {
  v <- tryCatch(isTRUE(expr), error = function(e) {
    message(id, ": ", conditionMessage(e)); FALSE })
  res[[id]] <<- v
}
n <- 10
# student
chk("distributions:4 integrate(dstudent_t)", isTRUE(all.equal(
  integrate(function(x) d_("student", x, list(mu = 10, sigma = 5,
                                              nu = 15)), -100, 100)$value, 1)))
chk("distributions:5 dstudent_t identity", isTRUE(all.equal(
  d_("student", 1, list(mu = 0, sigma = 5, nu = 10)), dt(1 / 5, 10) / 5)))
chk("distributions:9 rstudent_t length", length(r_("student", 10,
  list(mu = rnorm(10), sigma = 1:10, nu = 10))) == 10)
# von_mises
chk("distributions:50 dvon_mises length", length(d_("von_mises",
  runif(n, -pi, pi), list(mu = 1, kappa = 1:n))) == n)
chk("distributions:54 rvon_mises length", length(r_("von_mises", n,
  list(mu = rnorm(n), kappa = 0:(n - 1)))) == n)
# skew_normal
x <- rnorm(n, 10, 3)
chk("distributions:61 dskew_normal length", length(d_("skew_normal", x,
  list(mu = 1, sigma = 2, alpha = 1))) == n)
chk("distributions:73 rskew_normal length", length(r_("skew_normal", n,
  list(mu = rnorm(n), sigma = 10, alpha = -4:5))) == n)
# exgaussian
chk("distributions:80 dexgaussian length", length(d_("exgaussian", x,
  list(mu = 1, sigma = 2, beta = 1))) == n)
chk("distributions:85 rexgaussian length", length(r_("exgaussian", n,
  list(mu = rnorm(n), sigma = 10, beta = 1:10))) == n)
# inverse gaussian
xg <- rgamma(n, 10, 3)
chk("distributions:109 dinv_gaussian length", length(d_("inverse.gaussian",
  xg, list(mu = 1, shape = 1))) == n)
chk("distributions:113 rinv_gaussian length", length(r_("inverse.gaussian",
  n, list(mu = abs(rnorm(n)), shape = 1:10))) == n)
# beta_binomial
xp <- rpois(n, 1)
chk("distributions:123 dbeta_binomial length", length(d_("beta_binomial",
  pmin(xp, 2), list(mu = 0.4, phi = 1), list(trials = c(2, 10)))) == n)
chk("distributions:127 rbeta_binomial length", length(r_("beta_binomial",
  n, list(mu = 0.4, phi = 1), list(trials = c(2, 10)))) == n)
# asym_laplace
chk("distributions:146 dasym_laplace length", length(d_("asym_laplace", x,
  list(mu = 1, sigma = 2, quantile = 0.5))) == n)
chk("distributions:152 rasym_laplace length", length(r_("asym_laplace", n,
  list(mu = rnorm(n), sigma = 10, quantile = runif(n)))) == n)
# zero-inflated
chk("dzero_inflated_poisson length", length(d_("zero_inflated_poisson", xp,
  list(mu = 1, zi = 0.1))) == n)
chk("dzero_inflated_negbinomial length", length(d_(
  "zero_inflated_negbinomial", xp, list(mu = 2, shape = 5, zi = 0.1))) == n)
chk("dzero_inflated_binomial length", length(d_("zero_inflated_binomial",
  pmin(xp, 2), list(mu = 0.4, zi = 0.1), list(trials = c(2, 10)))) == n)
chk("dzero_inflated_beta_binomial length", length(d_(
  "zero_inflated_beta_binomial", pmin(xp, 2), list(mu = 0.4, phi = 1,
  zi = 0.1), list(trials = c(2, 10)))) == n)
xb <- c(rbeta(n - 2, 2, 3), 0, 0)
chk("dzero_inflated_beta length", length(d_("zero_inflated_beta", xb,
  list(mu = 2 / 5, phi = 5, zi = 0.1))) == n)
# hurdle
chk("dhurdle_poisson length", length(d_("hurdle_poisson", xp,
  list(mu = 1, hu = 0.1))) == n)
chk("dhurdle_negbinomial length", length(d_("hurdle_negbinomial", xp,
  list(mu = 2, shape = 5, hu = 0.1))) == n)
chk("dhurdle_gamma length", length(d_("hurdle_gamma", xp,
  list(mu = 3, shape = 1, hu = 0.1))) == n)
chk("dhurdle_lognormal length", length(d_("hurdle_lognormal", xp,
  list(mu = 2, sigma = 5, hu = 0.1))) == n)
# value checks against brms's own exported functions where defined
v <- c(
  zip = max(abs(d_("zero_inflated_poisson", xp, list(mu = 1, zi = 0.1)) -
    brms::dzero_inflated_poisson(xp, lambda = 1, zi = 0.1))),
  hup = max(abs(d_("hurdle_poisson", xp, list(mu = 1, hu = 0.1)) -
    brms::dhurdle_poisson(xp, lambda = 1, hu = 0.1))),
  hug = max(abs(d_("hurdle_gamma", xp, list(mu = 3, shape = 1, hu = 0.1)) -
    brms::dhurdle_gamma(xp, shape = 1, scale = 3, hu = 0.1))),
  sn = max(abs(d_("skew_normal", x, list(mu = 1, sigma = 2, alpha = 1)) -
    brms::dskew_normal(x, mu = 1, sigma = 2, alpha = 1))),
  vm = max(abs(d_("von_mises", xb, list(mu = 1, kappa = 2)) -
    brms::dvon_mises(xb, mu = 1, kappa = 2))))
r <- unlist(res)
cat(sprintf("%-44s %s\n", names(r), ifelse(r, "holds", "FAILS")), sep = "")
cat("\nholds:", sum(r), "of", length(r), "\n")
cat("max abs difference against brms's exported d function:\n")
print(signif(v, 3))
