# The family's own mean at the extreme nodes: hurdle at low mean, large sd.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(167)
d <- expand.grid(t = 1:6, g = factor(1:25)); d <- d[order(d$g, d$t), ]; n <- nrow(d)
d$x <- rnorm(n); e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
d$hp <- ifelse(rbinom(n, 1, .3) == 1, 0, pmax(rpois(n, exp(.8 + e)), 1))
fh <- frm(hp ~ x + ar(t, g, cov = TRUE), family = hurdle_poisson(), data = d)
gh <- fi("response_mean_gh")
fam <- fh$spec$responses[[1]]$family
for (s in c(1, 2.3, 4, 6)) for (e0 in c(-2, 1)) {
  v <- gh(e0, s^2, fam, list(mu = 0, hu = .3))
  x <- integrate(function(u) .7 * exp(e0 + s * u) / (-expm1(-exp(e0 + s * u))) * dnorm(u), -12, 12, rel.tol = 1e-10, subdivisions = 2000L)$value
  cat(sprintf("hurdle_poisson eta %+d sd %.1f: gh %.8g integrate %.8g rel %.2e\n", e0, s, v, x, v / x - 1))
}
