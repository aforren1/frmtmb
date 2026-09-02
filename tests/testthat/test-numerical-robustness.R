# Extreme-eta robustness of the family log-densities, and the huber()
# family.
#
# The at-the-optimum agreement tests cannot see any of this: a fit that
# converges never visits |eta| = 30. An optimizer step that overshoots
# does, and so does a separated predictor, a wide quadrature node and a
# frm_sample() tail draw. What follows tapes each log-density the way
# R/objective.R does - the dpar through its inverse link, the linear
# predictor alongside it under `.eta_<dpar>` - and asks for a finite
# value and a finite, correct gradient out there.

# The dpar list build_objective() hands a log-density, at the given
# per-dpar linear predictors.
dpars_at <- function(fam, etas) {
  out <- list()
  for (nm in names(etas)) {
    out[[nm]] <- fam$links[[nm]]$linkinv(etas[[nm]])
    out[[paste0(".eta_", nm)]] <- etas[[nm]]
  }
  out
}

# The taped negative log-density as a function of ONE dpar's linear
# predictor, every other dpar held at its own.
tape_nll <- function(fam, y, aterms, base_eta, sweep, extra = NULL) {
  f <- function(p) {
    etas <- base_eta
    etas[[sweep]] <- p[1] + 0 * base_eta[[sweep]]
    dp <- dpars_at(fam, etas)
    ll <- if (is.null(fam$extra_pars)) {
      fam$lpdf(y, dp, aterms)
    } else {
      fam$lpdf(y, dp, aterms, extra)
    }
    -sum(ll)
  }
  f
}

# Finite value, finite gradient, and a gradient that agrees with a
# central difference of the same function wherever that difference is
# itself computable.
expect_robust_at <- function(fam, y, aterms, base_eta, sweep, e0,
                             extra = NULL, tol = 1e-5) {
  f <- tape_nll(fam, y, aterms, base_eta, sweep, extra)
  tp <- RTMB::MakeTape(f, e0)
  v <- tp(e0)
  g <- as.numeric(tp$jacobian(e0))
  label <- paste0(fam$family, " ", sweep, " at eta = ", e0)
  expect_true(is.finite(v), label = paste("finite value:", label))
  expect_true(all(is.finite(g)), label = paste("finite gradient:", label))
  h <- 1e-4 * max(1, abs(e0))
  fd <- (f(e0 + h) - f(e0 - h)) / (2 * h)
  if (is.finite(fd)) {
    expect_lt(abs(g - fd) / max(1, abs(fd)), tol)
  }
  invisible(v)
}

ETA <- c(-30, 30)

test_that("binomial-family log-densities survive a separated predictor", {
  y01 <- c(0, 1)
  ybin <- c(0, 2, 5)
  for (e0 in ETA) {
    expect_robust_at(bernoulli(), y01, list(), list(mu = 0), "mu", e0)
    expect_robust_at(frmtmb:::fam_binomial(), ybin, list(trials = 5),
                     list(mu = 0), "mu", e0)
    expect_robust_at(beta_binomial(), ybin, list(trials = 5),
                     list(mu = 0, phi = log(5)), "mu", e0)
    expect_robust_at(Beta(), c(0.2, 0.5, 0.9), list(),
                     list(mu = 0, phi = log(5)), "mu", e0)
  }
})

test_that("cloglog saturates at single digits and the robust form does not", {
  # 1 - exp(-exp(eta)) is exactly 1 in double precision from eta = 4, so
  # this link's failure is not an exotic regime at all
  expect_equal(1 - exp(-exp(5)), 1)
  for (e0 in c(-30, -5, 5, 30)) {
    expect_robust_at(bernoulli("cloglog"), c(0, 1), list(), list(mu = 0),
                     "mu", e0)
    expect_robust_at(frmtmb:::fam_binomial("cloglog"), c(0, 2, 5),
                     list(trials = 5), list(mu = 0), "mu", e0)
  }
})

test_that("count families survive an underflowed mean", {
  ycnt <- c(0, 1, 3)
  for (e0 in ETA) {
    expect_robust_at(negbinomial(), ycnt, list(),
                     list(mu = 0, shape = log(2)), "mu", e0)
    expect_robust_at(nbinom1(), ycnt, list(), list(mu = 0, phi = log(2)),
                     "mu", e0)
    expect_robust_at(geometric(), ycnt, list(), list(mu = 0), "mu", e0)
  }
  # the two families whose negative-binomial SIZE does not move with the
  # mean stay finite even past the point where exp(eta) itself has
  # under- or overflowed. nbinom1's size is mu / phi, so it cannot: at
  # eta = -750 the size underflows to zero and the density degenerates,
  # which is a property of that parameterization rather than of the
  # arithmetic, and it behaved identically before this change.
  for (e0 in c(-750, 750)) {
    expect_robust_at(negbinomial(), ycnt, list(),
                     list(mu = 0, shape = log(2)), "mu", e0)
    expect_robust_at(geometric(), ycnt, list(), list(mu = 0), "mu", e0)
  }
})

test_that("zero-inflation and hurdle gates survive a separated gate", {
  ycnt <- c(0, 1, 3)
  for (e0 in ETA) {
    expect_robust_at(zero_inflated_poisson(), ycnt, list(),
                     list(mu = 0, zi = 0), "zi", e0)
    expect_robust_at(zero_inflated_negbinomial(), ycnt, list(),
                     list(mu = 0, shape = log(2), zi = 0), "zi", e0)
    expect_robust_at(zero_inflated_binomial(), c(0, 2, 5),
                     list(trials = 5), list(mu = 0, zi = 0), "zi", e0)
    expect_robust_at(zero_inflated_binomial(), c(0, 2, 5),
                     list(trials = 5), list(mu = 0, zi = 0), "mu", e0)
    expect_robust_at(zero_inflated_beta(), c(0, 0.4, 0.8), list(),
                     list(mu = 0, phi = log(5), zi = 0), "zi", e0)
    expect_robust_at(zero_inflated_asym_laplace(), c(0, 1, -1), list(),
                     list(mu = 0, sigma = 0, quantile = 0, zi = 0), "zi",
                     e0)
    expect_robust_at(hurdle_poisson(), ycnt, list(), list(mu = 0, hu = 0),
                     "hu", e0)
    expect_robust_at(hurdle_gamma(), c(0, 1, 2), list(),
                     list(mu = 0, shape = 0, hu = 0), "hu", e0)
    expect_robust_at(hurdle_lognormal(), c(0, 1, 2), list(),
                     list(mu = 0, sigma = 0, hu = 0), "hu", e0)
    # the zero-truncated Poisson normalizer log(1 - exp(-mu)) is the
    # other half of the hurdle family's exposure
    expect_robust_at(hurdle_poisson(), ycnt, list(), list(mu = 0, hu = 0),
                     "mu", e0)
    expect_robust_at(asym_laplace(), c(-1, 0, 1), list(),
                     list(mu = 0, sigma = 0, quantile = 0), "quantile", e0)
  }
})

test_that("ordinal log-densities survive extreme and near-coincident cuts", {
  yord <- c(1, 2, 3)
  wide <- list(tau_raw = c(-1, log(2)))
  close <- list(tau_raw = c(-1, log(1e-6)))
  for (e0 in ETA) {
    expect_robust_at(cumulative(), yord, list(), list(mu = 0), "mu", e0,
                     extra = wide)
    expect_robust_at(cumulative(), yord, list(), list(mu = 0), "mu", e0,
                     extra = close)
    expect_robust_at(sratio(), yord, list(), list(mu = 0), "mu", e0,
                     extra = wide)
    expect_robust_at(cratio(), yord, list(), list(mu = 0), "mu", e0,
                     extra = list(tau_raw = c(-1, 1)))
    expect_robust_at(acat(), yord, list(), list(mu = 0), "mu", e0,
                     extra = list(tau_raw = c(-1, 1)))
  }
})

test_that("multinomial-logit denominators accumulate in log space", {
  catf <- frmtmb:::fam_categorical_impl(c("mu2", "mu3"))
  for (e0 in c(ETA, 750)) {
    expect_robust_at(catf, c(1, 2, 3), list(), list(mu2 = 0, mu3 = 0),
                     "mu2", e0)
    expect_robust_at(acat(), c(1, 2, 3), list(), list(mu = 0), "mu", e0,
                     extra = list(tau_raw = c(-1, 1)))
  }
})

test_that("the robust forms leave a moderate-eta log-density where it was", {
  # the reference values are the SATURATING forms evaluated at ordinary
  # linear predictors, so a fit that used to converge somewhere still
  # converges there
  y01 <- c(0, 1)
  mu <- plogis(0.7)
  expect_equal(sum(RTMB::dbinom_robust(y01, 1, 0.7, log = TRUE)),
               sum(stats::dbinom(y01, 1, mu, log = TRUE)),
               tolerance = 1e-14)
  expect_equal(sum(RTMB::dnbinom_robust(c(0, 3), 0.7, 1.4 - log(2),
                                        log = TRUE)),
               sum(stats::dnbinom(c(0, 3), mu = exp(0.7), size = 2,
                                  log = TRUE)),
               tolerance = 1e-14)
  # log(p) / log(1 - p) off the log-odds, against the direct form
  for (e0 in c(-3, -0.5, 0, 0.5, 3)) {
    expect_equal(frmtmb:::log_inv_logit(e0), log(stats::plogis(e0)),
                 tolerance = 1e-14)
    expect_equal(frmtmb:::log1m_inv_logit(e0),
                 log1p(-stats::plogis(e0)), tolerance = 1e-13)
  }
  # cloglog's log-odds, against the round trip where the round trip is
  # still trustworthy
  for (e0 in c(-3, -0.5, 0, 1)) {
    p <- 1 - exp(-exp(e0))
    expect_equal(frmtmb:::frmtmb_links$cloglog$logit_eta(e0),
                 log(p / (1 - p)), tolerance = 1e-10)
  }
})

test_that("the eta-scale pass-through does not change a fit", {
  set.seed(303)
  n <- 150
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(1:15, each = 10)))
  d$yb <- stats::rbinom(n, 1, stats::plogis(0.4 + 0.9 * d$x))
  fit <- frm(bf(yb ~ x + (1 | g)), family = bernoulli(), data = d)
  ref <- lme4::glmer(yb ~ x + (1 | g), d, family = stats::binomial())
  expect_lt(max(abs(fixef(fit)$mu - lme4::fixef(ref))), 1e-3)
  expect_lt(abs(as.numeric(stats::logLik(fit)) -
                  as.numeric(stats::logLik(ref))), 1e-4)
})

test_that("separation is fitted instead of crashing", {
  # perfectly separated data: the ML estimate runs off to infinity, so
  # the optimizer WILL walk into the saturated region. It has to come
  # back with numbers.
  d <- data.frame(x = c(-3, -2, -1, 1, 2, 3), y = c(0, 0, 0, 1, 1, 1))
  fit <- suppressWarnings(frm(bf(y ~ x), family = bernoulli(), data = d))
  expect_true(all(is.finite(fixef(fit)$mu)))
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # the fitted slope is large, and the log-likelihood is essentially 0
  expect_gt(fixef(fit)$mu[["x"]], 5)
  expect_lt(abs(as.numeric(stats::logLik(fit))), 1e-3)
})

# --- huber() ---------------------------------------------------------

test_that("huber()'s normalizer is the density's own", {
  for (k in c(0.5, 1.345, 2, 5)) {
    num <- stats::integrate(function(u) exp(-frmtmb:::huber_rho(u, k)),
                            -Inf, Inf, rel.tol = 1e-12)$value
    expect_equal(frmtmb:::huber_norm(k), num, tolerance = 1e-10)
    # and the density integrates to one
    mass <- stats::integrate(
      function(u) exp(-frmtmb:::huber_rho(u, k)) / frmtmb:::huber_norm(k),
      -Inf, Inf, rel.tol = 1e-12)$value
    expect_equal(mass, 1, tolerance = 1e-10)
    # the closed-form second moment, for pearson residuals
    m2 <- stats::integrate(
      function(u) u^2 * exp(-frmtmb:::huber_rho(u, k)) /
        frmtmb:::huber_norm(k), -Inf, Inf, rel.tol = 1e-12)$value
    expect_equal(frmtmb:::huber_var_u(k), m2, tolerance = 1e-10)
  }
})

test_that("the branch-free rho is the piecewise one", {
  k <- 1.345
  u <- c(-6, -1.5, -k, -0.4, 0, 0.4, k, 1.5, 6)
  ref <- ifelse(abs(u) <= k, u^2 / 2, k * abs(u) - k^2 / 2)
  expect_equal(frmtmb:::huber_rho(u, k), ref, tolerance = 1e-15)
})

test_that("huber() estimates match Huber's own estimating equations", {
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1.5 + 0.8 * d$x - 0.4 * d$z + stats::rnorm(n)
  d$y[1:10] <- d$y[1:10] + 25          # gross contamination
  k <- 1.345
  fit <- suppressWarnings(frm(bf(y ~ x + z), family = huber(), data = d))
  X <- stats::model.matrix(~ x + z, d)
  u <- (d$y - as.vector(X %*% fixef(fit)$mu)) / sigma(fit)
  psi <- pmin(pmax(u, -k), k)
  # the score of the Huber log-likelihood, derived by hand: rho' is the
  # clamped residual, and d/d log sigma gives sum(u psi(u)) = n
  expect_lt(max(abs(crossprod(X, psi))), 1e-3)
  expect_lt(abs(sum(u * psi) - n), 1e-3)

  # and against a hand-rolled optimization of the same log-likelihood
  lz <- log(frmtmb:::huber_norm(k))
  nll <- function(p) {
    s <- exp(p[4])
    r <- (d$y - as.vector(X %*% p[1:3])) / s
    sum(log(s) + lz + ifelse(abs(r) <= k, r^2 / 2, k * abs(r) - k^2 / 2))
  }
  ref <- stats::nlminb(c(stats::coef(stats::lm(y ~ x + z, d)),
                         log(stats::mad(d$y))), nll,
                       control = list(rel.tol = 1e-15, x.tol = 1e-14,
                                      iter.max = 5000, eval.max = 20000))
  # rho's kink stops any optimizer short of machine precision, so the
  # sharp statement is that our objective is no worse than the
  # reference's, not that the coordinates agree to 1e-8
  expect_lt(-as.numeric(stats::logLik(fit)) - ref$objective, 1e-8)
  expect_lt(max(abs(fixef(fit)$mu - ref$par[1:3])), 1e-4)
  expect_lt(abs(sigma(fit) - exp(ref$par[4])), 1e-4)
})

test_that("huber() tracks MASS::rlm, and the gap is the scale", {
  skip_if_not_installed("MASS")
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1.5 + 0.8 * d$x - 0.4 * d$z + stats::rnorm(n)
  rl <- MASS::rlm(y ~ x + z, data = d, psi = MASS::psi.huber, k = 1.345)
  fit <- suppressWarnings(frm(bf(y ~ x + z), family = huber(), data = d))
  # rlm holds the scale at a MAD-type estimate and iterates the
  # location; huber() estimates sigma by ML jointly with mu, so the
  # coefficients agree only to about 1e-2
  expect_lt(max(abs(fixef(fit)$mu - stats::coef(rl))), 2e-2)
  # hold sigma where rlm holds it and the two estimators coincide,
  # which is what says the difference is entirely the scale
  fixed <- suppressWarnings(
    frm(bf(y ~ x + z, sigma = rl$s), family = huber(), data = d))
  expect_lt(max(abs(fixef(fixed)$mu - stats::coef(rl))), 1e-4)
})

test_that("huber() collapses to gaussian() as k grows", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(150))
  d$y <- 2 - 0.6 * d$x + stats::rnorm(150)
  g <- frm(bf(y ~ x), family = gaussian(), data = d)
  h <- suppressWarnings(frm(bf(y ~ x), family = huber(k = 20), data = d))
  expect_lt(max(abs(fixef(h)$mu - fixef(g)$mu)), 1e-4)
  expect_lt(abs(sigma(h) - sigma(g)), 1e-4)
  expect_lt(abs(as.numeric(stats::logLik(h)) -
                  as.numeric(stats::logLik(g))), 1e-6)
})

test_that("huber() bounds the influence of an outlier", {
  set.seed(12)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- 1 + 0.8 * d$x + stats::rnorm(120)
  clean <- suppressWarnings(frm(bf(y ~ x), family = huber(), data = d))
  d2 <- d
  d2$y[1:4] <- d2$y[1:4] + 40
  dirty <- suppressWarnings(frm(bf(y ~ x), family = huber(), data = d2))
  gdirty <- frm(bf(y ~ x), family = gaussian(), data = d2)
  gclean <- frm(bf(y ~ x), family = gaussian(), data = d)
  # the gaussian slope moves several times as far as the huber one
  expect_lt(abs(fixef(dirty)$mu[["x"]] - fixef(clean)$mu[["x"]]),
            abs(fixef(gdirty)$mu[["x"]] - fixef(gclean)$mu[["x"]]))
})

test_that("huber()'s simulator draws from huber()'s density", {
  set.seed(9)
  for (k in c(0.8, 3)) {
    u <- frmtmb:::rhuber_u(2e5, k)
    expect_lt(abs(mean(u)), 0.02)
    expect_lt(abs(stats::var(u) / frmtmb:::huber_var_u(k) - 1), 0.02)
    inside <- stats::integrate(
      function(v) exp(-frmtmb:::huber_rho(v, k)) / frmtmb:::huber_norm(k),
      -k, k)$value
    expect_lt(abs(mean(abs(u) < k) - inside), 0.01)
  }
})

test_that("huber() plumbs through the fit methods", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(80))
  d$y <- 2 - 0.6 * d$x + stats::rnorm(80)
  fit <- suppressWarnings(frm(bf(y ~ x), family = huber(), data = d))
  expect_equal(unname(fitted(fit)),
               unname(predict(fit, type = "response")))
  expect_true(all(is.finite(residuals(fit, type = "pearson"))))
  expect_true(all(is.finite(residuals(fit, type = "deviance"))))
  expect_length(simulate(fit)[[1]], 80)
  expect_identical(family(fit)$family, "huber")
  # k is a constant of the family, not a dpar
  expect_identical(family(fit)$dpars, c("mu", "sigma"))
  expect_error(huber(k = 0), "tuning constant")
  expect_error(huber(k = c(1, 2)), "tuning constant")
  # reachable by name, as every registry family is
  expect_identical(frmtmb:::as_frmtmb_family("huber")$family, "huber")
})
