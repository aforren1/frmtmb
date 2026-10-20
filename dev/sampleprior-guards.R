# Lane sampleprior, fix round 1, item 5: the support guards added to
# prior_base_logdens() change nothing inside the support. In one process
# (dev/rtmb-pitfalls.md, item 21), each comparison runs with the guarded
# densities and then with 0.71.0's unguarded ones put back into the
# namespace: MAP fits with each bounded density, the tape's value and
# gradient at points inside the support, and a seeded frm_sample() run
# under the defaults. Data seed 20261023.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "\n")
core <- asNamespace("frmtmb")
guarded <- core$prior_base_logdens
# 0.71.0's densities, the guards taken out again
unguarded <- function(x, dist) {
  switch(dist$kind,
    exponential = log(dist$rate) - dist$rate * x,
    gamma = dist$shape * log(dist$rate) - lgamma(dist$shape) +
      (dist$shape - 1) * log(x) - dist$rate * x,
    inv_gamma = dist$shape * log(dist$scale) - lgamma(dist$shape) -
      (dist$shape + 1) * log(x) - dist$scale / x,
    beta = (dist$shape1 - 1) * log(x) + (dist$shape2 - 1) * log1p(-x) -
      (lgamma(dist$shape1) + lgamma(dist$shape2) -
         lgamma(dist$shape1 + dist$shape2)),
    lognormal = RTMB::dnorm(log(x), dist$location, dist$scale,
                            log = TRUE) - log(x),
    guarded(x, dist))
}
with_dens <- function(f, expr) {
  old <- core$prior_base_logdens
  utils::assignInNamespace("prior_base_logdens", f, "frmtmb")
  on.exit(utils::assignInNamespace("prior_base_logdens", old, "frmtmb"))
  expr
}
set.seed(20261023)
n <- 80
x <- stats::rnorm(n)
dnb <- data.frame(x = x, y = stats::rnbinom(n, mu = exp(1 + 0.3 * x),
                                            size = 2))
dzi <- data.frame(x = x, y = ifelse(stats::runif(n) < 0.3, 0,
                                    stats::rpois(n, exp(1 + 0.3 * x))))
dg <- data.frame(x = x, y = 1 + x + stats::rnorm(n))
cases <- list(
  list("negbinomial, shape ~ gamma(1, 0.5)", bf(y ~ x), negbinomial(), dnb,
       set_prior("gamma(1, 0.5)", class = "shape")),
  list("negbinomial, shape ~ inv_gamma(0.4, 0.3)", bf(y ~ x), negbinomial(),
       dnb, set_prior("inv_gamma(0.4, 0.3)", class = "shape")),
  list("negbinomial, shape ~ lognormal(0, 1)", bf(y ~ x), negbinomial(),
       dnb, set_prior("lognormal(0, 1)", class = "shape")),
  list("zip, zi ~ beta(1, 1)", bf(y ~ x), zero_inflated_poisson(), dzi,
       set_prior("beta(1, 1)", class = "zi")),
  list("zip, zi ~ beta(2, 5)", bf(y ~ x), zero_inflated_poisson(), dzi,
       set_prior("beta(2, 5)", class = "zi")),
  list("gaussian, sigma ~ exponential(1)", bf(y ~ x), gaussian(), dg,
       set_prior("exponential(1)", class = "sigma")))
for (cs in cases) {
  f1 <- suppressWarnings(frm(cs[[2]], family = cs[[3]], data = cs[[4]],
                             prior = cs[[5]]))
  f0 <- with_dens(unguarded, suppressWarnings(frm(cs[[2]], family = cs[[3]],
                                                  data = cs[[4]],
                                                  prior = cs[[5]])))
  e1 <- unlist(f1$estimates)
  e0 <- unlist(f0$estimates)
  # the tape of each, at 20 points around the estimate
  p <- f1$obj$par
  set.seed(1)
  pts <- replicate(20, p + stats::rnorm(length(p), 0, 0.3))
  dv <- dg_ <- 0
  for (k in seq_len(ncol(pts))) {
    v1 <- f1$obj$fn(pts[, k])
    v0 <- f0$obj$fn(pts[, k])
    g1 <- f1$obj$gr(pts[, k])
    g0 <- f0$obj$gr(pts[, k])
    dv <- max(dv, abs(v1 - v0) / abs(v0))
    dg_ <- max(dg_, max(abs(g1 - g0)) / max(abs(g0)))
  }
  cat(sprintf(paste0("%-40s MAP estimates identical: %s (max rel diff ",
                     "%.2g) | tape at 20 points: max rel diff value %.2g, ",
                     "gradient %.2g\n"), cs[[1]], identical(e1, e0),
              max(abs(e1 - e0) / pmax(abs(e0), 1e-300)), dv, dg_))
}
# a seeded run under the defaults, whose zi takes beta(1, 1) and shape
# inv_gamma(0.4, 0.3): the guarded and the unguarded objectives
for (cs in list(list("zinb defaults", bf(y ~ x),
                     zero_inflated_negbinomial(), dzi))) {
  run <- function() {
    as.matrix(suppressWarnings(suppressMessages(frm_sample(
      cs[[2]], family = cs[[3]], data = cs[[4]], chains = 2, iter = 1000,
      refresh = 0, seed = 7))))
  }
  m1 <- run()
  m0 <- with_dens(unguarded, run())
  cat(sprintf(paste0("%-40s draws identical: %s; max |diff| %.3g; ",
                     "posterior means %s against %s\n"), cs[[1]],
              identical(m1, m0), max(abs(m1 - m0)),
              paste(sprintf("%.4f", colMeans(m1)[1:4]), collapse = " "),
              paste(sprintf("%.4f", colMeans(m0)[1:4]), collapse = " ")))
}
# the control: the swap reaches the tape. beta(1, 1) on a coefficient
# (with brms's lb and ub left out on purpose, which warns), at a value
# outside (0, 1): the guarded tape is NaN, 0.71.0's is finite
pr <- suppressWarnings(set_prior("beta(1, 1)", class = "b", coef = "x"))
f1 <- suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = dg,
                           prior = pr))
f0 <- with_dens(unguarded, suppressWarnings(frm(bf(y ~ x),
                                                family = gaussian(),
                                                data = dg, prior = pr)))
p <- f1$obj$par
p[2] <- 1.5
cat(sprintf("control, beta(1, 1) on b_x at 1.5: guarded %s, unguarded %s\n",
            format(f1$obj$fn(p)), format(f0$obj$fn(p))))
# fix round 2: the logit link's Jacobian in natural_log_jacobian(),
# against 0.71.0's log(mu_eta) put back, on the zero-inflation fits:
# equal to rounding where p (1 - p) does not underflow
old_jac <- function(x, link) log(abs(link$mu_eta(x)))
for (cs in cases[4:5]) {
  f1 <- suppressWarnings(frm(cs[[2]], family = cs[[3]], data = cs[[4]],
                             prior = cs[[5]]))
  oj <- core$natural_log_jacobian
  utils::assignInNamespace("natural_log_jacobian", old_jac, "frmtmb")
  f0 <- suppressWarnings(frm(cs[[2]], family = cs[[3]], data = cs[[4]],
                             prior = cs[[5]]))
  utils::assignInNamespace("natural_log_jacobian", oj, "frmtmb")
  e1 <- unlist(f1$estimates)
  e0 <- unlist(f0$estimates)
  cat(sprintf("%-40s logit Jacobian: MAP max rel diff %.2g, nll %.17g vs %.17g\n",
              cs[[1]], max(abs(e1 - e0) / pmax(abs(e0), 1e-300)),
              f1$opt$objective, f0$opt$objective))
}
