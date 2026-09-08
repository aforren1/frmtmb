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

# --- brms link parity -------------------------------------------------
#
# probit, probit_approx, cauchit, softit, softplus, squareplus, sqrt,
# log1p and `1/mu^2` complete the registry against brms 2.23.0. The
# saturation points asserted below are measured, not assumed; each one
# is the linear predictor at which the named path first passes 1e-8
# relative error against an independent reference.

NEW_LINKS <- c("probit", "probit_approx", "cauchit", "softit",
               "softplus", "squareplus", "sqrt", "log1p", "1/mu^2")

# an eta each link is meaningful on: sqrt and `1/mu^2` need eta > 0 for
# a positive mean, the rest take the whole line
link_etas <- function(nm) {
  switch(nm, sqrt = c(0.4, 1.1, 2.3), `1/mu^2` = c(0.3, 1.1, 4),
         log1p = c(-0.7, 0.3, 1.2), c(-2, -0.5, 0.3, 1.7))
}

test_that("the new links reproduce brms's own definitions", {
  skip_if_not_installed("brms")
  for (nm in NEW_LINKS) {
    lk <- frmtmb:::frmtmb_links[[nm]]
    e <- link_etas(nm)
    # brms's R-side inv_link() answers pnorm() for probit_approx while
    # its Stan program uses Phi_approx(); the Stan form is the one a
    # ported model was fitted with, so that is what is asserted.
    want <- if (identical(nm, "probit_approx")) {
      1 / (1 + exp(-(0.07056 * e^3 + 1.5976 * e)))
    } else {
      brms:::inv_link(e, nm)
    }
    expect_equal(lk$linkinv(e), want, tolerance = 1e-14,
                 label = paste0(nm, " linkinv"))
    mu <- if (identical(nm, "probit_approx")) stats::pnorm(e) else want
    expect_equal(lk$linkfun(mu), brms:::link(mu, nm), tolerance = 1e-14,
                 label = paste0(nm, " linkfun"))
  }
})

test_that("every registry field tapes with a finite value and gradient", {
  # a link field that cannot be taped is a link that cannot be fitted,
  # and the robust fields are read from inside the log-density
  start <- list(logit = 0.63, probit = 0.63, probit_approx = 0.63,
                cauchit = 0.63, cloglog = 0.63, softit = 0.63,
                logm1 = 2.4, power12 = 1.4, log1p = 0.4, identity = 0.4)
  for (nm in names(frmtmb:::frmtmb_links)) {
    lk <- frmtmb:::frmtmb_links[[nm]]
    e0 <- link_etas(nm)[1]
    if (nm %in% c("inverse", "log", "sqrt", "1/mu^2")) e0 <- 0.7
    for (f in c("linkfun", "linkinv", "mu_eta", "logit_eta", "log_eta")) {
      if (is.null(lk[[f]])) next
      x0 <- if (identical(f, "linkfun")) start[[nm]] %||% 0.9 else e0
      tp <- RTMB::MakeTape(function(p) sum(lk[[f]](p)), x0)
      expect_true(is.finite(tp(x0)),
                  label = paste(nm, f, "taped value"))
      expect_true(all(is.finite(as.numeric(tp$jacobian(x0)))),
                  label = paste(nm, f, "taped gradient"))
    }
  }
})

test_that("mu_eta is the derivative predict(se.fit) thinks it is", {
  skip_if_not_installed("numDeriv")
  for (nm in NEW_LINKS) {
    lk <- frmtmb:::frmtmb_links[[nm]]
    e <- link_etas(nm)
    nd <- vapply(e, function(z) numDeriv::grad(lk$linkinv, z), 0)
    expect_equal(lk$mu_eta(e), nd, tolerance = 1e-8,
                 label = paste0(nm, " mu_eta"))
  }
  # RTMB::pnorm's own AD derivative is the normal density exactly, which
  # is what probit's mu_eta is written as
  for (e0 in c(-3, -0.5, 0, 1.2, 4)) {
    tp <- RTMB::MakeTape(function(p) sum(RTMB::pnorm(p)), e0)
    expect_equal(as.numeric(tp$jacobian(e0)),
                 frmtmb:::frmtmb_links$probit$mu_eta(e0),
                 tolerance = 1e-15)
  }
})

test_that("probit's log-odds holds where the round trip has saturated", {
  probit <- frmtmb:::frmtmb_links$probit
  # agreement with the round trip while the round trip is still
  # trustworthy
  for (e0 in c(-3, -0.5, 0, 0.5, 3)) {
    p <- stats::pnorm(e0)
    expect_equal(probit$logit_eta(e0), log(p / (1 - p)), tolerance = 1e-12)
  }
  # the plain path dies at eta = 6.4: 1 - pnorm(eta) rounds to zero once
  # pnorm(eta) rounds to one
  expect_equal(1 - stats::pnorm(9), 0)
  expect_identical(log1p(-stats::pnorm(9)), -Inf)
  # the robust field is still exact there, and stays so until pnorm
  # itself underflows near |eta| = 38.2
  for (e0 in c(7, 9, 20, 37)) {
    expect_equal(probit$logit_eta(e0),
                 stats::pnorm(e0, log.p = TRUE) -
                   stats::pnorm(-e0, log.p = TRUE),
                 tolerance = 1e-12)
    expect_true(is.finite(probit$logit_eta(e0)))
    expect_true(is.finite(probit$logit_eta(-e0)))
  }
  expect_false(is.finite(probit$logit_eta(40)))
})

test_that("probit_approx's log-odds is the cubic and never saturates", {
  pa <- frmtmb:::frmtmb_links$probit_approx
  for (e0 in c(-3, -0.5, 0, 0.5, 3)) {
    p <- pa$linkinv(e0)
    expect_equal(pa$logit_eta(e0), log(p / (1 - p)), tolerance = 1e-12)
  }
  # the plain path is past 1e-8 by eta = 5.9 and exactly saturated by 9
  expect_identical(log1p(-pa$linkinv(9)), -Inf)
  # the inverse link is a logistic OF the cubic, so the log-odds is the
  # cubic at any eta at all
  for (e0 in c(9, 100, 1e5)) {
    expect_equal(pa$logit_eta(e0), 0.07056 * e0^3 + 1.5976 * e0,
                 tolerance = 1e-14)
    expect_true(is.finite(pa$logit_eta(e0)))
  }
})

test_that("softit's log-odds is the softplus's log", {
  s <- frmtmb:::frmtmb_links$softit
  for (e0 in c(-3, -0.5, 0, 0.5, 3)) {
    p <- s$linkinv(e0)
    expect_equal(s$logit_eta(e0), log(p / (1 - p)), tolerance = 1e-12)
    expect_equal(s$logit_eta(e0), log(log1p(exp(e0))), tolerance = 1e-12)
  }
  # softit's upper tail is polynomial, so the plain path survives much
  # further than the logit's: mu rounds to one only past eta = 1e16,
  # where the robust field is still exact
  expect_identical(log1p(-s$linkinv(1e17)), -Inf)
  expect_true(is.finite(frmtmb:::log1m_inv_logit(s$logit_eta(1e17))))
})

test_that("cauchit needs no robust field and says so by having none", {
  cc <- frmtmb:::frmtmb_links$cauchit
  expect_null(cc$logit_eta)
  expect_null(cc$log_eta)
  # the Cauchy tail is polynomial: 1 - mu is 1 / (pi eta) to leading
  # order, so the plain round trip is still good to 1e-10 relative at a
  # linear predictor no optimizer will ever reach
  for (e0 in c(1e3, 1e6)) {
    expect_equal(log1p(-cc$linkinv(e0)),
                 stats::pcauchy(e0, lower.tail = FALSE, log.p = TRUE),
                 tolerance = 1e-10)
  }
})

test_that("squareplus's log mean is asinh, where the sum has cancelled", {
  sq <- frmtmb:::frmtmb_links$squareplus
  for (e0 in c(-8, -1, 0, 1, 8)) {
    expect_equal(sq$log_eta(e0), log(sq$linkinv(e0)), tolerance = 1e-13)
  }
  # (eta + sqrt(eta^2 + 4)) / 2 cancels for eta below about -1e5 and is
  # exactly zero by -1e9, where the true log mean is an ordinary -20.7
  expect_equal(sq$linkinv(-1e9), 0)
  expect_identical(log(sq$linkinv(-1e9)), -Inf)
  expect_equal(sq$log_eta(-1e9), -log(1e9), tolerance = 1e-9)
  expect_true(is.finite(sq$log_eta(-1e12)))
})

test_that("softplus and sqrt log_eta pick the robust density branch", {
  # Neither field beats the plain round trip at computing log(mu). What
  # each does is switch the family onto its robust log-density: without
  # a log_eta, negbinomial forms `mu + mu^2 / shape`, whose excess over
  # mu is lost to rounding once mu is tiny, and the log-density then
  # carries noise a central difference can see.
  sp <- frmtmb:::frmtmb_links$softplus
  expect_equal(sp$log_eta(-700), log(sp$linkinv(-700)), tolerance = 1e-12)
  expect_equal(sp$log_eta(-700), -700, tolerance = 1e-12)
  expect_identical(sp$log_eta(-746), -Inf)
  for (e0 in c(-3, 0, 2, 30)) {
    expect_equal(sp$log_eta(e0), log(sp$linkinv(e0)), tolerance = 1e-13)
  }
  sq <- frmtmb:::frmtmb_links$sqrt
  for (e0 in c(-30, -1, 0.5, 30)) {
    expect_equal(sq$log_eta(e0), log(sq$linkinv(e0)), tolerance = 1e-13)
  }
  # log(eta^2), not 2 * log(abs(eta)): the same number with no kink on
  # the tape at zero
  tp <- RTMB::MakeTape(function(p) sum(sq$log_eta(p)), 0.5)
  expect_equal(as.numeric(tp$jacobian(0.5)), 4, tolerance = 1e-12)
})

test_that("the links with no exact robust form leave the field absent", {
  # `1/mu^2` never loses a digit through the round trip, and
  # inverse.gaussian has no robust branch to select either way
  expect_null(frmtmb:::frmtmb_links[["1/mu^2"]]$log_eta)
  expect_equal(log(frmtmb:::frmtmb_links[["1/mu^2"]]$linkinv(1e300)),
               -0.5 * log(1e300), tolerance = 1e-14)
  # log1p sits on a dpar that may be negative, so neither field applies
  expect_null(frmtmb:::frmtmb_links$log1p$log_eta)
  expect_null(frmtmb:::frmtmb_links$log1p$logit_eta)
})

test_that("the new (0, 1) links survive a separated predictor", {
  y01 <- c(0, 1)
  ybin <- c(0, 2, 5)
  for (lk in c("probit", "probit_approx", "cauchit", "softit")) {
    for (e0 in ETA) {
      # bernoulli and binomial take the log-odds straight into
      # dbinom_robust and never exponentiate it back, so any magnitude
      # of log-odds is fine here
      expect_robust_at(bernoulli(lk), y01, list(), list(mu = 0), "mu", e0)
      expect_robust_at(frmtmb:::fam_binomial(lk), ybin, list(trials = 5),
                       list(mu = 0), "mu", e0)
      # beta's shapes are mu * phi and (1 - mu) * phi, so
      # dpar_complement() has to exponentiate the log-odds back and a
      # shape underflows to zero past |log-odds| = 745. probit reaches
      # only 454 at eta = 30, but probit_approx's log-odds is CUBIC in
      # eta and is already at 1953 there, so it is swept where the
      # shape is still representable.
      be <- if (identical(lk, "probit_approx")) sign(e0) * 20 else e0
      expect_robust_at(Beta(lk), c(0.2, 0.5, 0.9), list(),
                       list(mu = 0, phi = log(5)), "mu", be)
    }
  }
  # and that limit is a property of the beta shapes, not of the link
  expect_equal(exp(frmtmb:::log_inv_logit(
    frmtmb:::frmtmb_links$probit_approx$logit_eta(-30))), 0)
})

test_that("the new positive-mean links survive an extreme predictor", {
  ycnt <- c(0, 1, 3)
  for (lk in c("softplus", "squareplus", "sqrt")) {
    for (e0 in c(-30, 30)) {
      expect_robust_at(negbinomial(lk), ycnt, list(),
                       list(mu = 0, shape = log(2)), "mu", e0)
    }
  }
})

test_that("the new links fit a GLM and match stats::glm", {
  set.seed(404)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$nt <- sample(4:10, n, TRUE)
  # both optimizers tightened: at their defaults the two stop at points
  # 2e-6 apart in the coefficients while agreeing to 1e-9 in the
  # log-likelihood, which measures the stopping rules and not the links
  ctl <- frmtmb_control(optCtrl = list(rel.tol = 1e-14, x.tol = 1e-12,
                                       eval.max = 5000, iter.max = 5000))
  gctl <- stats::glm.control(epsilon = 1e-14, maxit = 200)
  for (lk in c("probit", "cauchit")) {
    p <- if (identical(lk, "probit")) stats::pnorm(0.3 + 0.8 * d$x)
         else stats::pcauchy(0.3 + 0.8 * d$x)
    d$y <- stats::rbinom(n, d$nt, p)
    fit <- suppressWarnings(frm(bf(y | trials(nt) ~ x + z),
                                family = stats::binomial(link = lk),
                                data = d, control = ctl))
    g <- stats::glm(cbind(y, nt - y) ~ x + z, data = d,
                    family = stats::binomial(link = lk), control = gctl)
    expect_lt(max(abs(fixef(fit)$mu - stats::coef(g))), 1e-6)
    expect_lt(abs(as.numeric(stats::logLik(fit)) -
                    as.numeric(stats::logLik(g))), 1e-6)
  }
  # drawn through the link, with sqrt(mu) bounded away from zero: IRLS
  # on a sqrt link diverges as soon as a step puts the linear predictor
  # negative, which a mean near zero invites
  d$yp <- stats::rpois(n, (2 + 0.4 * pmax(pmin(d$x, 3), -3))^2)
  fp <- suppressWarnings(frm(bf(yp ~ x),
                             family = stats::poisson(link = "sqrt"),
                             data = d, control = ctl))
  gp <- stats::glm(yp ~ x, data = d, start = fixef(fp)$mu,
                   family = stats::poisson(link = "sqrt"), control = gctl)
  expect_lt(max(abs(fixef(fp)$mu - stats::coef(gp))), 1e-6)

  # `1/mu^2` is the inverse Gaussian canonical link and stats spells it
  # the same way brms does, so this one comparison is available. The
  # mean is drawn THROUGH the link, with a linear predictor bounded well
  # away from zero: IRLS on this link takes sqrt(eta) and diverges the
  # moment a step puts eta below it.
  d$eta_i <- 2 + 0.3 * pmax(pmin(d$x, 3.5), -3.5)
  d$yi <- stats::rgamma(n, 8, 8 * sqrt(d$eta_i))
  fi <- suppressWarnings(frm(bf(yi ~ x),
                             family = stats::inverse.gaussian(),
                             data = d, control = ctl))
  expect_identical(family(fi)$links$mu$name, "1/mu^2")
  gi <- suppressWarnings(stats::glm(
    yi ~ x, data = d, family = stats::inverse.gaussian(),
    start = fixef(fi)$mu, control = gctl))
  expect_lt(max(abs(fixef(fi)$mu - stats::coef(gi))), 1e-6)
})

test_that("predict(se.fit) is the delta method through every new link", {
  skip_if_not_installed("numDeriv")
  set.seed(77)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * d$x))
  nd <- data.frame(x = c(-1.5, -0.4, 0.6, 2))
  for (lk in c("probit", "probit_approx", "cauchit", "softit")) {
    fit <- suppressWarnings(frm(bf(y ~ x), family = bernoulli(lk),
                                data = d))
    pl <- stats::predict(fit, newdata = nd, type = "link", se.fit = TRUE)
    pr <- stats::predict(fit, newdata = nd, type = "response",
                         se.fit = TRUE)
    lo <- frmtmb:::frmtmb_links[[lk]]
    dn <- vapply(pl$fit, function(z) numDeriv::grad(lo$linkinv, z), 0)
    expect_equal(unname(pr$se.fit), unname(abs(dn) * pl$se.fit),
                 tolerance = 1e-7, label = paste(lk, "se.fit"))
    # and the bands conditional_effects draws stay inside the support
    ce <- conditional_effects(fit)[[1]]
    expect_true(all(is.finite(ce$estimate__)))
    expect_true(all(ce$lower__ > 0 & ce$upper__ < 1))
  }
})

test_that("the links stats::make.link rejects still reach frmtmb", {
  # stats::binomial() and stats::Gamma() validate their link string
  # through make.link() before frmtmb ever sees it, so a brms link
  # stats has never heard of has to come in through one of frmtmb's own
  # family constructors.
  expect_error(stats::binomial(link = "softit"), "not recognised")
  set.seed(5)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * d$x))
  for (lk in c("probit_approx", "softit")) {
    fit <- suppressWarnings(frm(bf(y ~ x), family = bernoulli(lk),
                                data = d))
    expect_identical(family(fit)$links$mu$name, lk)
    expect_true(all(is.finite(fixef(fit)$mu)))
  }
  d$yw <- stats::rweibull(n, 2, exp(0.5 + 0.3 * d$x))
  for (lk in c("softplus", "squareplus")) {
    fit <- suppressWarnings(frm(bf(yw ~ x), family = weibull(lk), data = d))
    expect_identical(family(fit)$links$mu$name, lk)
    expect_true(all(is.finite(fixef(fit)$mu)))
  }
})

# --- conditional_effects() bands through every new link ---------------
#
# The band is built on the scale it is symmetric on and then moved to
# the response scale. Two things go wrong there and both are asserted
# here rather than assumed: an inverse link that DECREASES sends the
# lower linear predictor to the upper response bound, and a link whose
# domain does not reach an end has no bound at all on that side.
# `1/mu^2` is where this lane first met both, because
# stats::inverse.gaussian()'s default link only became fittable when the
# registry gained the name.

# one family per new link that accepts it, and the data drawn through
# that link so the fit sits where the link is well conditioned
band_case <- function(lk, seed = 5, n = 200) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  e <- 0.3 + 0.5 * d$x
  fam <- switch(lk,
    probit = , probit_approx = , cauchit = , softit = {
      lo <- frmtmb:::frmtmb_links[[lk]]
      d$y <- stats::rbinom(n, 1, lo$linkinv(e))
      bernoulli(lk)
    },
    softplus = , squareplus = , sqrt = , log1p = {
      lo <- frmtmb:::frmtmb_links[[lk]]
      # a strictly positive mean: sqrt and log1p need the predictor kept
      # off their turning points
      ep <- 1.5 + 0.4 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rpois(n, lo$linkinv(ep))
      negbinomial(lk)
    },
    `1/mu^2` = {
      ei <- 2 + 0.3 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rgamma(n, 8, 8 * sqrt(ei))
      stats::inverse.gaussian()
    },
    inverse = {
      ei <- 2 + 0.3 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rgamma(n, 8, 8 * ei)
      stats::Gamma(link = "inverse")
    },
    stop("no case for ", lk))
  list(fit = suppressWarnings(frm(bf(y ~ x), family = fam, data = d)),
       data = d)
}

test_that("conditional_effects bands are ordered through every new link", {
  for (lk in c("probit", "probit_approx", "cauchit", "softit",
               "softplus", "squareplus", "sqrt", "log1p", "1/mu^2")) {
    cs <- band_case(lk)
    expect_identical(family(cs$fit)$links$mu$name, lk,
                     label = paste(lk, "link reached the fit"))
    r <- suppressWarnings(conditional_effects(cs$fit)[[1]])
    ok <- is.finite(r$lower__) & is.finite(r$upper__)
    # NaN is never a bound: a bound the link cannot reach is NA
    expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)),
                 label = paste(lk, "band has no NaN"))
    expect_true(all(r$lower__[ok] <= r$upper__[ok]),
                label = paste(lk, "band is ordered"))
    # and the estimate sits inside its own band
    expect_true(all(r$estimate__[ok] >= r$lower__[ok] &
                      r$estimate__[ok] <= r$upper__[ok]),
                label = paste(lk, "estimate inside the band"))
    expect_true(any(ok), label = paste(lk, "band is not all NA"))
  }
})

test_that("a decreasing link does not invert the band", {
  # Gamma(inverse) inverted on 100 of 100 rows before ce_band_ends():
  # 1/eta decreases, so linkinv sent the lower linear predictor to the
  # upper response bound and the two were never swapped back. This is
  # the pre-existing half of the defect - it reproduces with the link
  # registry untouched - and `1/mu^2` is the half that was new.
  for (lk in c("inverse", "1/mu^2")) {
    cs <- band_case(lk)
    r <- suppressWarnings(conditional_effects(cs$fit)[[1]])
    ok <- is.finite(r$lower__) & is.finite(r$upper__)
    expect_gt(sum(ok), 0)
    expect_equal(sum(r$lower__[ok] > r$upper__[ok]), 0,
                 label = paste(lk, "inverted rows"))
    expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)))
  }
})

test_that("a bound outside the link's domain is NA, once, and named", {
  # 1/sqrt(eta) has no value at eta <= 0, so a band wide enough to reach
  # zero has no UPPER response bound - the mean runs away as eta falls
  # to zero. The bound must be NA on that side, not NaN, and not filled
  # in from the other end.
  set.seed(11)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  mig <- 1 / sqrt(pmax(0.5 + 0.2 * d$x, 0.05))
  d$y <- stats::rgamma(n, shape = 20, rate = 20 / mig)
  fit <- suppressWarnings(frm(bf(y ~ x),
                              family = stats::inverse.gaussian(),
                              data = d))
  expect_warning(r <- conditional_effects(fit)[[1]],
                 "1/mu\\^2.*does not reach the band")
  ok <- is.finite(r$lower__) & is.finite(r$upper__)
  expect_equal(sum(r$lower__[ok] > r$upper__[ok]), 0)
  expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)))
  # the missing side is the UPPER one, and the lower bound survives
  gone <- !is.finite(r$upper__)
  expect_gt(sum(gone), 0)
  expect_true(all(is.finite(r$lower__[gone])))
  expect_true(all(is.na(r$upper__[gone])))
})

test_that("ce_band_ends orders, refuses a pole and keeps the sides", {
  # the helper on its own, where the cases can be built exactly
  inc <- frmtmb:::frmtmb_links$log
  dec <- frmtmb:::frmtmb_links[["1/mu^2"]]
  # increasing: ends arrive in order and stay
  b <- frmtmb:::ce_band_ends(inc$linkinv, c(-1, 0), c(1, 2),
                             inc$mu_eta(c(0, 1)))
  expect_equal(b$lower, exp(c(-1, 0)))
  expect_equal(b$upper, exp(c(1, 2)))
  expect_identical(b$outside, 0L)
  # decreasing: the ends are swapped, not sorted
  b2 <- frmtmb:::ce_band_ends(dec$linkinv, c(1, 4), c(4, 9),
                              dec$mu_eta(c(2, 6)))
  expect_equal(b2$lower, 1 / sqrt(c(4, 9)))
  expect_equal(b2$upper, 1 / sqrt(c(1, 4)))
  expect_identical(b2$outside, 0L)
  # an end the link cannot map: NA on the side it belongs to
  b3 <- frmtmb:::ce_band_ends(dec$linkinv, c(-1), c(4), dec$mu_eta(2))
  expect_equal(b3$lower, 0.5)
  expect_true(is.na(b3$upper))
  expect_identical(b3$outside, 1L)
  # a pole BETWEEN two finite ends: neither bounds the response, so
  # both are refused. 1/eta from -1 to 1 passes through infinity.
  invl <- frmtmb:::frmtmb_links$inverse
  b4 <- frmtmb:::ce_band_ends(invl$linkinv, -1, 1, invl$mu_eta(-1))
  expect_true(is.na(b4$lower))
  expect_true(is.na(b4$upper))
  expect_identical(b4$outside, 1L)
  # mixed-sign derivative in ONE call: sqrt's mu_eta is 2 * eta, so the
  # link decreases left of zero and increases right of it. The swap is
  # per row, so both halves come back ordered; deciding it once for the
  # vector would invert whichever half lost the vote.
  sq <- frmtmb:::frmtmb_links$sqrt
  lo <- c(-3, 1); hi <- c(-1, 3)
  b5 <- frmtmb:::ce_band_ends(sq$linkinv, lo, hi, sq$mu_eta((lo + hi) / 2))
  expect_equal(b5$lower, c(1, 1))
  expect_equal(b5$upper, c(9, 9))
  expect_true(all(b5$lower <= b5$upper))
  expect_identical(b5$outside, 0L)
  # and a length-1 derivative against a longer band swaps every row
  b6 <- frmtmb:::ce_band_ends(dec$linkinv, c(1, 4), c(4, 9), dec$mu_eta(2))
  expect_equal(b6$lower, 1 / sqrt(c(4, 9)))
  expect_equal(b6$upper, 1 / sqrt(c(1, 4)))
})

test_that("predict(se.fit) reports a standard error, not an interval", {
  # H1 is a band defect and predict() has no band: it returns
  # |mu_eta| * se_eta, which is non-negative by construction and has no
  # endpoints to order. Asserted so the claim is measured, not assumed.
  for (lk in c("1/mu^2", "inverse")) {
    cs <- band_case(lk)
    nd <- data.frame(x = c(-1.5, -0.4, 0.6, 2))
    p <- stats::predict(cs$fit, newdata = nd, type = "response",
                        se.fit = TRUE)
    expect_named(p, c("fit", "se.fit"))
    expect_true(all(is.finite(p$se.fit)))
    expect_true(all(p$se.fit >= 0))
  }
})
