# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived, and this file is a main contributor
# to CRAN-condition check time.
skip_on_cran()

# trunc() must reach the post-fit surface, not just the likelihood:
# fitted(), predict(type = "response"), residuals() and simulate() all
# describe the TRUNCATED distribution. [brms#1923, #1903]

# Brute-force E[Y | lb <= Y <= ub] over the support.
trunc_mean_disc <- function(dfun, lb, ub) {
  ys <- seq(lb, min(ub, 5000))
  p <- dfun(ys)
  sum(ys * p) / sum(p)
}
trunc_mean_cont <- function(dfun, lb, ub) {
  num <- stats::integrate(function(y) y * dfun(y), lb, ub,
                          rel.tol = 1e-11)$value
  den <- stats::integrate(dfun, lb, ub, rel.tol = 1e-11)$value
  num / den
}

fit_trunc_pois <- function() {
  set.seed(1)
  y <- stats::rpois(6000, 3)
  y <- y[y >= 2 & y <= 6][1:400]
  frm(bf(y | trunc(lb = 2, ub = 6) ~ 1) + poisson(),
      data = data.frame(y = y))
}

test_that("truncated poisson fitted() is the truncated mean", {
  fit <- fit_trunc_pois()
  mu <- exp(fixef(fit)$mu[[1]])
  want <- trunc_mean_disc(function(y) stats::dpois(y, mu), 2, 6)
  expect_equal(unname(fitted(fit)[1]), want, tolerance = 1e-10)
  expect_equal(unname(predict(fit, type = "response")[1]), want,
               tolerance = 1e-10)
  # the fitted mean is well away from the untruncated mu: this is the
  # whole defect (2.976 vs 3.375)
  expect_gt(abs(want - mu), 0.3)
  # dpar-scale predictions describe the LATENT parameter and stay
  # untruncated
  expect_equal(unname(predict(fit, type = "conditional")[1]), mu,
               tolerance = 1e-10)
  expect_equal(unname(predict(fit, dpar = "mu", type = "response")[1]), mu,
               tolerance = 1e-10)
  expect_equal(unname(exp(predict(fit, type = "link")[1])), mu,
               tolerance = 1e-10)
  # an intercept-only correctly specified fit reproduces the sample mean
  expect_equal(mean(residuals(fit)), 0, tolerance = 1e-8)
})

test_that("truncated gaussian and lognormal means match integration", {
  set.seed(11)
  n <- 1500
  x <- stats::rnorm(n)
  ystar <- 1 + 0.5 * x + stats::rnorm(n)
  keep <- ystar > 0.5
  gd <- data.frame(y = ystar[keep], x = x[keep])
  gfit <- frm(bf(y | trunc(lb = 0.5) ~ x) + gaussian(), data = gd)
  dp <- frmtmb:::eval_dpars(gfit)[["y"]]
  ref <- vapply(seq_len(nrow(gd)), function(i) {
    m <- dp$mu[i]; s <- dp$sigma[i]
    trunc_mean_cont(function(z) stats::dnorm(z, m, s), 0.5, m + 15 * s)
  }, 0)
  expect_vector_equal(as.numeric(fitted(gfit)), ref, tol = 1e-6)
  # the truncation-aware residuals of a correctly specified model centre
  expect_lt(abs(mean(residuals(gfit))), 0.02)

  set.seed(12)
  ln <- stats::rlnorm(4000, 0.4, 0.8)
  ln <- ln[ln >= 1 & ln <= 5][1:300]
  lfit <- frm(bf(y | trunc(lb = 1, ub = 5) ~ 1) + lognormal(),
              data = data.frame(y = ln))
  ld <- frmtmb:::eval_dpars(lfit)[["y"]]
  lref <- trunc_mean_cont(
    function(z) stats::dlnorm(z, ld$mu[1], ld$sigma[1]), 1, 5)
  expect_equal(unname(fitted(lfit)[1]), lref, tolerance = 1e-8)
  expect_lt(abs(mean(residuals(lfit))), 0.02)
})

test_that("simulate() respects trunc() bounds and hits the right mean", {
  fit <- fit_trunc_pois()
  sims <- as.matrix(simulate(fit, nsim = 60, seed = 7))
  expect_true(all(sims >= 2 & sims <= 6))
  expect_equal(mean(sims), unname(fitted(fit)[1]), tolerance = 0.05)

  # bounds that exclude essentially all the fitted mass cannot be filled
  # by rejection; the error must say so rather than spin
  fam <- frmtmb:::fam_poisson()
  expect_error(
    frmtmb:::sim_response(fam, list(mu = rep(3, 20)),
                          list(trunc_lb = rep(400, 20)), 20,
                          max_iter = 3L),
    "rejection sampling"
  )
})

test_that("DHARMa residuals on a truncated fit are uniform", {
  skip_if_not_installed("DHARMa")
  fit <- fit_trunc_pois()
  dh <- dharma_residuals(fit, nsim = 250, seed = 5)
  expect_gt(stats::ks.test(dh$scaledResiduals, "punif")$p.value, 0.01)
})

test_that("OSA residuals on truncated fits are standard normal", {
  fit <- fit_trunc_pois()
  expect_gt(stats::ks.test(residuals(fit, type = "osa"),
                           "pnorm")$p.value, 0.01)
  set.seed(11)
  n <- 1500
  x <- stats::rnorm(n)
  ystar <- 1 + 0.5 * x + stats::rnorm(n)
  keep <- ystar > 0.5
  gd <- data.frame(y = ystar[keep], x = x[keep])
  gfit <- frm(bf(y | trunc(lb = 0.5) ~ x) + gaussian(), data = gd)
  osa <- residuals(gfit, type = "osa")
  expect_gt(stats::ks.test(osa, "pnorm")$p.value, 0.01)
  # ... and equal the analytic truncated-normal PIT
  dp <- frmtmb:::eval_dpars(gfit)[["y"]]
  Fl <- stats::pnorm((0.5 - dp$mu) / dp$sigma)
  ref <- stats::qnorm(
    (stats::pnorm((gd$y - dp$mu) / dp$sigma) - Fl) / (1 - Fl))
  expect_vector_equal(osa, ref, tol = 1e-4)
})

test_that("truncation bounds follow newdata", {
  set.seed(21)
  x <- stats::rnorm(800)
  y <- stats::rpois(800, exp(0.9 + 0.4 * x))
  keep <- y >= 1 & y <= 8
  dd <- data.frame(y = y[keep], x = x[keep])
  # literal bounds carry over to any newdata
  fit <- frm(bf(y | trunc(lb = 1, ub = 8) ~ x) + poisson(), data = dd)
  nd <- data.frame(x = c(-1, 0, 1))
  mu <- exp(predict(fit, newdata = nd, type = "link"))
  want <- vapply(mu, function(m) {
    trunc_mean_disc(function(z) stats::dpois(z, m), 1, 8)
  }, 0)
  expect_vector_equal(predict(fit, newdata = nd, type = "response"), want,
                      tol = 1e-10)
  # in-sample newdata reproduces fitted()
  expect_vector_equal(predict(fit, newdata = dd, type = "response"),
                      as.numeric(fitted(fit)), tol = 1e-8)

  # a bound given as a VARIABLE is data, so newdata must supply it
  dd$lo <- 1
  dd$hi <- 8
  vfit <- frm(bf(y | trunc(lb = lo, ub = hi) ~ x) + poisson(), data = dd)
  expect_vector_equal(as.numeric(fitted(vfit)), as.numeric(fitted(fit)),
                      tol = 1e-8)
  nd2 <- data.frame(x = c(-1, 0, 1), lo = 1, hi = 8)
  expect_vector_equal(predict(vfit, newdata = nd2, type = "response"),
                      want, tol = 1e-10)
  expect_error(predict(vfit, newdata = nd, type = "response"),
               "trunc\\(lb = lo\\)")
})

test_that("zero-truncated poisson matches glmmTMB likelihood and mean", {
  skip_if_not_installed("glmmTMB")
  set.seed(31)
  x <- stats::rnorm(500)
  y <- stats::rpois(500, exp(0.7 + 0.4 * x))
  keep <- y > 0
  dd <- data.frame(y = y[keep], x = x[keep])
  fit <- frm(bf(y | trunc(lb = 1) ~ x) + poisson(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, data = dd,
                          family = glmmTMB::truncated_poisson())
  expect_loglik_equal(fit, ref, tol = 1e-5)
  # glmmTMB's response prediction IS the truncated mean, so it doubles as
  # a reference for ours; brute force pins both
  mu <- exp(as.numeric(predict(fit, type = "link")))
  brute <- vapply(mu, function(m) {
    trunc_mean_disc(function(z) stats::dpois(z, m), 1, Inf)
  }, 0)
  ours <- as.numeric(predict(fit, type = "response"))
  expect_vector_equal(ours, brute, tol = 1e-8)
  expect_vector_equal(ours, unname(stats::predict(ref, type = "response")),
                      tol = 1e-4)
  # mu/(1 - exp(-mu)) is the closed form; it is not mu
  expect_gt(max(abs(ours - mu)), 0.1)
})
