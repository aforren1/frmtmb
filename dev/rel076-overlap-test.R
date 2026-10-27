## ---- the 0.76.0 consolidation: the round's three changes together -----
# mo()'s coding from the minimum, the residual-correlation priors (cosy
# sampled on brms's [0, 1] through one logistic, sigma's prior read at
# the innovation sd under a first-order cov = TRUE ar()) and log_lik()
# of a residual covariance term were each tested without the others.
# These blocks sample real chains and read their draws: the mo() mean is
# written from brms's definition, D * sum(simo[1:(k - min)]), and the
# reference is the dense conditional of each group at each draw, from
# the draws' own columns (dev/rel076-overlap.R is the probe).

test_that("cosy() and mo() together: log_lik() of the draws is the dense
           conditional, with brms's mo() coding and cosy on [0, 1]", {
  skip_sampler()
  set.seed(20261027)
  d <- expand.grid(t = 1:4, g = 1:30)
  d$x <- stats::rnorm(nrow(d))
  d$k <- sample(1:4, nrow(d), replace = TRUE)
  d$y <- 1 + 0.5 * d$x + c(0, 0.6, 0.9, 1.5)[d$k] +
    stats::rnorm(30, 0, 0.4)[d$g] + stats::rnorm(nrow(d), 0, 0.9)
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ x + mo(k) + cosy(t, g), data = d, chains = 1, iter = 400,
    seed = 3, refresh = 0)), c("Effective Samples Size", "R-hat"))
  m <- as.matrix(ds$draws)
  sx <- grep("^simo_mok1\\[", colnames(m))
  # values 1 to 4: three weights, as brms codes them
  expect_length(sx, 3L)
  expect_true(all(m[, "cosy"] >= 0 & m[, "cosy"] <= 1))
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(nrow(m), nrow(d)))
  for (s in seq(1, nrow(m), by = 20)) {
    cs <- c(0, cumsum(m[s, sx]))
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x +
      m[s, "bsp_mok"] * length(sx) * cs[d$k - min(d$k) + 1]
    rho <- m[s, "cosy"]
    ref <- lla_grouped(d, mu, rep(m[s, "sigma"], nrow(d)),
                       function(L) ifelse(L == 0, 1, rho))
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
  skip_if_not_installed("loo")
  lo <- suppressWarnings(loo(ds))
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
  expect_length(lo$pointwise[, "elpd_loo"], nrow(d))
})

test_that("ar(cov = TRUE) under the default priors: log_lik() of the
           draws reads sigma as the marginal sd", {
  skip_sampler()
  # the sigma prior is read at the innovation sd, sigma sqrt(1 - ar^2),
  # but the draws' sigma column is frmtmb's marginal sd, and the
  # conditional is of the marginal covariance sigma^2 ar^|lag|
  set.seed(20261027)
  d <- expand.grid(t = 1:6, g = 1:8)
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(1:8, function(g) {
    as.numeric(stats::arima.sim(list(ar = 0.9), 6, sd = 0.45))
  }))
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ x + ar(t, g, cov = TRUE), data = d, chains = 1, iter = 400,
    seed = 3, refresh = 0)), c("Effective Samples Size", "R-hat"))
  m <- as.matrix(ds$draws)
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(nrow(m), nrow(d)))
  for (s in seq(1, nrow(m), by = 20)) {
    a1 <- m[s, "ar[1]"]
    ref <- lla_grouped(d, m[s, "b_Intercept"] + m[s, "b_x"] * d$x,
                       rep(m[s, "sigma"], nrow(d)), function(L) a1^L)
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
  skip_if_not_installed("loo")
  lo <- suppressWarnings(loo(ds))
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
})
