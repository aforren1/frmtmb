# conditional_effects(method = "predict") on a maximum-likelihood fit
# reports estimate__ as the MEDIAN of the plug-in predictive draws, the
# same draws its band is a quantile band of, as brms reports the median
# of its predictive draws (posterior_summary(robust = TRUE)). Up to
# 0.67.0 it reported the expected response, which differs on every
# skewed or discrete family (dev/fixes-ce-pred.R: Poisson at x = -1, 0,
# 1, brms 1 2 4 against frmtmb 1.46 2.52 4.35). BREAKING at 0.68.0, the
# user's decision of 2026-10-06.

cpm_ce <- function(fit, ndraws, ...) {
  set.seed(11)
  conditional_effects(fit, effects = "x", method = "predict",
                      resolution = 5, ndraws = ndraws, ...)$x
}

# the Monte Carlo sd of a sample median of `n` draws from a density
# `f0` at its median: 1 / (2 f0 sqrt(n))
cpm_mc_sd <- function(f0, n) 1 / (2 * f0 * sqrt(n))

test_that("a continuous family's estimate__ is the predictive median", {
  set.seed(1)
  d <- data.frame(x = rnorm(300))
  d$y <- exp(0.3 + 0.5 * d$x + rnorm(300))
  fit <- frm(y ~ x, family = lognormal(), data = d)
  nd <- 2001L
  ce <- cpm_ce(fit, nd)
  mu <- as.numeric(frm_linpred(fit, newdata = data.frame(x = ce$x),
                               dpar = "mu", type = "link"))
  s <- sigma(fit)
  med <- exp(mu)
  mean_y <- exp(mu + s^2 / 2)
  # within five Monte Carlo sds of the median, and far from the mean,
  # which is exp(s^2 / 2) = 1.6 times it here
  sd_med <- cpm_mc_sd(stats::dlnorm(med, mu, s), nd)
  expect_true(all(abs(ce$estimate__ - med) < 5 * sd_med))
  expect_true(all(abs(ce$estimate__ - mean_y) > 20 * sd_med))
  # the 50% point of the band's own distribution
  expect_true(all(ce$lower__ <= ce$estimate__ &
                    ce$estimate__ <= ce$upper__))
})

test_that("a discrete family's estimate__ is a whole count", {
  set.seed(2)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(0.8 + 0.6 * d$x))
  fit <- frm(y ~ x, family = poisson(), data = d)
  ce <- cpm_ce(fit, 401L)
  lam <- exp(as.numeric(frm_linpred(fit, newdata = data.frame(x = ce$x),
                                    dpar = "mu", type = "link")))
  expect_identical(ce$estimate__, round(ce$estimate__))
  # the exact median of a Poisson, up to one count where the CDF sits
  # near one half
  expect_true(all(abs(ce$estimate__ - stats::qpois(0.5, lam)) <= 1))
  expect_true(all(ce$lower__ <= ce$estimate__ &
                    ce$estimate__ <= ce$upper__))
  # binomial with trials() pinned: a count out of the pinned trials
  d$nt <- sample(5:20, 300, TRUE)
  d$k <- rbinom(300, d$nt, stats::plogis(0.3 + 0.5 * d$x))
  fb <- frm(k | trials(nt) ~ x, family = binomial(), data = d)
  ceb <- suppressMessages(cpm_ce(fb, 401L, conditions = list(nt = 10)))
  expect_identical(ceb$estimate__, round(ceb$estimate__))
  expect_true(all(ceb$estimate__ >= 0 & ceb$estimate__ <= 10))
})

test_that("a truncated family's estimate__ is the truncated median", {
  set.seed(4)
  dt <- data.frame(x = rnorm(400))
  dt$y <- 1 + 0.5 * dt$x + rnorm(400)
  dt <- dt[dt$y > 0.5, ]
  ft <- frm(y | trunc(lb = 0.5) ~ x, family = gaussian(), data = dt)
  nd <- 2001L
  ce <- cpm_ce(ft, nd)
  mu <- as.numeric(frm_linpred(ft, newdata = data.frame(x = ce$x),
                               dpar = "mu", type = "link"))
  s <- sigma(ft)
  p0 <- stats::pnorm(0.5, mu, s)
  med <- stats::qnorm(p0 + 0.5 * (1 - p0), mu, s)
  f0 <- stats::dnorm(med, mu, s) / (1 - p0)
  expect_true(all(ce$estimate__ >= 0.5))
  expect_true(all(abs(ce$estimate__ - med) < 5 * cpm_mc_sd(f0, nd)))
})

test_that("the epred display stays the expected response", {
  # only method = "predict" changed; frmtmb.sample's draws route reported
  # the median of its draws already and is unchanged
  set.seed(5)
  d <- data.frame(x = rnorm(200))
  d$y <- rpois(200, exp(0.5 + 0.4 * d$x))
  fit <- frm(y ~ x, family = poisson(), data = d)
  ce <- conditional_effects(fit, effects = "x", resolution = 5)$x
  ep <- as.numeric(frm_linpred(fit, newdata = data.frame(x = ce$x),
                               type = "response"))
  expect_lt(max(abs(ce$estimate__ - ep)), 1e3 * .Machine$double.eps *
              max(ep))
})

test_that("the median is brms's at the same parameters", {
  # brms's own conditional_effects(method = "posterior_predict") on a
  # brmsfit whose 4001 draws all sit at frmtmb's estimates is the same
  # plug-in predictive distribution; both report a Monte Carlo median
  # of it (dev/fixes-ce-pred.R measured the same with 4000 draws)
  skip_unless_brms_fit()
  set.seed(1)
  d <- data.frame(x = rnorm(300))
  d$y <- exp(0.3 + 0.5 * d$x + rnorm(300))
  d$k <- rpois(300, exp(0.8 + 0.6 * d$x))
  xs <- c(-1, 0, 1)
  nd <- 4001L
  for (case in c("lognormal", "poisson")) {
    f <- if (case == "poisson") k ~ x else y ~ x
    fam <- if (case == "poisson") poisson() else lognormal()
    fit <- frm(f, family = fam, data = d)
    bb <- brms_fixed_fit(brms::bf(f), get(case, asNamespace("brms"))(),
                         d, fit, ndraws = nd)
    set.seed(3)
    a <- suppressWarnings(brms::conditional_effects(
      bb, "x", method = "posterior_predict",
      int_conditions = list(x = xs)))[[1]]
    set.seed(3)
    b <- conditional_effects(fit, "x", method = "predict", ndraws = nd,
                             int_conditions = list(x = xs))[[1]]
    if (case == "poisson") {
      lam <- exp(as.numeric(frm_linpred(fit, newdata = data.frame(x = xs),
                                        dpar = "mu", type = "link")))
      # both are the exact median, a count, unless the CDF sits near one
      # half there; at these three points it does not
      expect_identical(b$estimate__, stats::qpois(0.5, lam))
      expect_identical(a$estimate__, b$estimate__)
    } else {
      mu <- as.numeric(frm_linpred(fit, newdata = data.frame(x = xs),
                                   dpar = "mu", type = "link"))
      sdm <- cpm_mc_sd(stats::dlnorm(exp(mu), mu, sigma(fit)), nd)
      # two independent Monte Carlo medians: five sds of their difference
      expect_true(all(abs(a$estimate__ - b$estimate__) < 5 * sqrt(2) * sdm))
    }
  }
})
