# xbeta(), zero_inflated_beta_binomial() and hurdle_cumulative(), brms
# 2.23.0's three families that frmtmb lacked, and cse(), brms's second
# name for cs().
#
# Every reference density here is written from stats::dbeta(),
# stats::pbeta(), lbeta(), lchoose() and the link CDFs after brms's Stan
# functions (chunks/fun_xbeta.stan, fun_zero_inflated_beta_binomial.stan
# and brms:::stan_hurdle_ordinal_lpmf()), not from the package code.
# Tolerances are ratios: to the double epsilon for a density, to the
# measured log-likelihood for a fit, and to the measured standard error
# for a coefficient. dev/fams2-validate.R prints the same comparisons.

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

ref_zibb <- function(y, size, mu, phi, zi) {
  a <- mu * phi
  b <- (1 - mu) * phi
  base <- lchoose(size, y) + lbeta(y + a, size - y + b) - lbeta(a, b)
  ifelse(y == 0, log(zi + (1 - zi) * exp(base)), log1p(-zi) + base)
}

# each link's distribution function and its survivor, both from the
# stats functions' own tails, so that a category probability deep in
# either tail is formed without cancelling
hc_cdf <- list(
  logit = list(F = stats::plogis,
               S = function(x) stats::plogis(x, lower.tail = FALSE)),
  probit = list(F = stats::pnorm,
                S = function(x) stats::pnorm(x, lower.tail = FALSE)),
  cauchit = list(F = stats::pcauchy,
                 S = function(x) stats::pcauchy(x, lower.tail = FALSE)),
  cloglog = list(F = function(x) -expm1(-exp(x)),
                 S = function(x) exp(-exp(x))))

ref_hc <- function(y, eta, hu, disc, tau, link) {
  cd <- hc_cdf[[link]]
  K1 <- length(tau)
  n <- max(length(y), length(eta))
  y <- rep(y, length.out = n)
  a <- disc * (c(tau, Inf)[pmin(pmax(y, 1), K1 + 1)] - eta)
  b <- disc * (c(-Inf, tau)[pmin(pmax(y, 1), K1 + 1)] - eta)
  # F(a) - F(b) from whichever tail both ends sit in
  p <- ifelse(a + b > 0, cd$S(b) - cd$S(a), cd$F(a) - cd$F(b))
  ifelse(y == 0, log(hu), log1p(-hu) + log(p))
}

rel_diff <- function(a, b) max(abs(a - b) / abs(b))
ULPS <- 64 * .Machine$double.eps

sim_xbeta <- function(seed, n = 500, ngrp = 0) {
  set.seed(seed)
  x <- rnorm(n)
  g <- factor(if (ngrp > 0) rep(seq_len(ngrp), length.out = n) else 1)
  u <- if (ngrp > 0) rnorm(ngrp, 0, 0.4)[g] else 0
  mu <- stats::plogis(0.3 + 0.6 * x + u)
  kap <- exp(-2 + 0.5 * x)
  z <- stats::rbeta(n, mu * 6, (1 - mu) * 6)
  data.frame(y = pmin(pmax((1 + 2 * kap) * z - kap, 0), 1), x = x, g = g)
}

sim_zibb <- function(seed, n = 500, ngrp = 0) {
  set.seed(seed)
  x <- rnorm(n)
  g <- factor(if (ngrp > 0) rep(seq_len(ngrp), length.out = n) else 1)
  u <- if (ngrp > 0) rnorm(ngrp, 0, 0.4)[g] else 0
  tr <- sample(4:20, n, TRUE)
  mu <- stats::plogis(-0.4 + 0.5 * x + u)
  yb <- stats::rbinom(n, tr, stats::rbeta(n, mu * 4, (1 - mu) * 4))
  data.frame(y = ifelse(runif(n) < stats::plogis(-1 + 0.4 * x), 0L, yb),
             x = x, tr = tr, g = g)
}

sim_hc <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  u <- stats::rlogis(n) + 0.8 * x
  y <- ifelse(runif(n) < stats::plogis(-0.6 + 0.5 * z), 0L,
              1L + (u > -1) + (u > 0.3) + (u > 1.5))
  data.frame(y = y, x = x, z = z)
}

# ------------------------------------------------------------ constructors

test_that("the constructors carry brms's names, dpars and default links", {
  f <- xbeta()
  expect_identical(f$family, "xbeta")
  expect_identical(f$dpars, c("mu", "phi", "kappa"))
  expect_identical(c(f$link, f$link_phi, f$link_kappa),
                   c("logit", "log", "log"))
  expect_identical(xbeta(cloglog)$link, "cloglog")
  g <- zero_inflated_beta_binomial()
  expect_identical(g$dpars, c("mu", "phi", "zi"))
  expect_identical(c(g$link, g$link_phi, g$link_zi),
                   c("logit", "log", "logit"))
  expect_identical(brmsfamily("zi_beta_binomial")$family,
                   "zero_inflated_beta_binomial")
  h <- hurdle_cumulative()
  expect_identical(h$type, "ordinal")
  expect_identical(h$dpars, c("mu", "hu", "disc"))
  expect_identical(c(h$link, h$link_hu, h$link_disc),
                   c("logit", "logit", "log"))
  expect_identical(hurdle_cumulative("cauchit")$link, "cauchit")
  expect_identical(brmsfamily("hu_cumulative", "probit")$link, "probit")
})

test_that("links, responses and options brms refuses are refused by name", {
  expect_error(xbeta("1/mu"), "not a supported link for family 'xbeta'",
               fixed = TRUE)
  expect_error(xbeta(link_phi = "sqrt"), "'sqrt' is not a supported link",
               fixed = TRUE)
  expect_error(zero_inflated_beta_binomial(link_zi = "log"),
               "not a supported link for parameter 'zi'", fixed = TRUE)
  expect_error(zero_inflated_beta_binomial("sqrt"),
               "not a supported link for family", fixed = TRUE)
  expect_error(hurdle_cumulative(link = "log"),
               "not a supported link for family 'hurdle_cumulative'",
               fixed = TRUE)
  expect_error(hurdle_cumulative(link_hu = "probit"),
               "not a supported link for parameter 'hu'", fixed = TRUE)
  expect_error(hurdle_cumulative(link_disc = "logit"),
               "not a supported link for parameter 'disc'", fixed = TRUE)
  expect_error(hurdle_cumulative(threshold = "equidistant"),
               "is not implemented", fixed = TRUE)
  expect_error(hurdle_cumulative(threshold = "wobbly"),
               "takes one of", fixed = TRUE)
  expect_error(frm(bf(y ~ x), family = xbeta(),
                   data = data.frame(y = c(0, 0.5, 1.2), x = 1:3)),
               "xbeta: response must be in [0, 1], where", fixed = TRUE)
  d <- data.frame(y = c(0, 3, 13), tr = c(10, 10, 10), x = 1:3)
  expect_error(frm(bf(y | trials(tr) ~ x),
                   family = zero_inflated_beta_binomial(), data = d),
               "response must be integer counts in [0, trials]",
               fixed = TRUE)
  expect_error(frm(bf(y ~ x), family = zero_inflated_beta_binomial(),
                   data = d), "trials")
  d <- data.frame(y = c(0, 1, -1, 2), x = 1:4)
  expect_error(frm(bf(y ~ x), family = hurdle_cumulative(), data = d),
               "requires either non-negative integers or ordered factors")
  d <- data.frame(y = c(0, 1, 0, 1), x = 1:4)
  expect_error(frm(bf(y ~ x), family = hurdle_cumulative(), data = d),
               "Could not extract the number of thresholds")
})

test_that("cens(), trunc(), cs(), thres(gr = ) and osa are refused", {
  d <- sim_zibb(1, n = 60)
  expect_error(frm(bf(y | trials(tr) + cens(x > 1) ~ x),
                   family = zero_inflated_beta_binomial(), data = d),
               "need a family with a CDF")
  dx <- sim_xbeta(2, n = 80)
  expect_error(frm(bf(y | trunc(ub = 0.9) ~ x), family = xbeta(), data = dx),
               "need a family with a CDF")
  dh <- sim_hc(3, n = 120)
  dh$g <- gl(2, 60)
  expect_error(frm(bf(y ~ cs(x)), family = hurdle_cumulative(), data = dh),
               "cs() needs an sratio, cratio, or acat family", fixed = TRUE)
  expect_error(frm(bf(y | thres(gr = g) ~ x), family = hurdle_cumulative(),
                   data = dh), "not thres(gr = )", fixed = TRUE)
  fit <- frm(bf(y ~ x), family = hurdle_cumulative(), data = dh)
  expect_error(residuals(fit, type = "osa"), "point mass", fixed = TRUE)
  fx <- frm(bf(y ~ x), family = xbeta(), data = dx)
  expect_error(residuals(fx, type = "osa"), "point mass", fixed = TRUE)
})

# ----------------------------------------------------------------- density

test_that("xbeta's density is brms's, at both ends and inside, on the tape", {
  fam <- xbeta()
  y <- c(0, 1, 1e-6, 0.3, 0.97)
  for (eta in c(-6, 0, 5)) for (lphi in c(-1, 2, 5)) for (lk in c(-7, 0)) {
    mu <- stats::plogis(eta)
    ref <- ref_xbeta(y, mu, exp(lphi), exp(lk))
    dp <- list(mu = mu, phi = exp(lphi), kappa = exp(lk))
    got <- fam$lpdf(y, lapply(dp, rep, length.out = 5L), list())
    # the ends go through log_ibeta_half(), which holds 5e-13 of
    # stats::pbeta() at these shapes (dev/fams2-p1-sweep-pkg.txt)
    expect_lt(rel_diff(got, ref), 1e3 * ULPS)
    F <- RTMB::MakeTape(function(p) {
      e <- p[1] + 0 * y
      sum(fam$lpdf(y, list(mu = 1 / (1 + exp(-e)), .eta_mu = e,
                           phi = exp(p[2]), kappa = exp(p[3])), list()))
    }, c(0, 0, 0))
    expect_lt(rel_diff(F(c(eta, lphi, lk)), sum(ref)), 1e3 * ULPS)
  }
})

test_that("kappa running to 0 warns, and nothing else does", {
  # PIN: without this check the fit ran kappa to exp(-22.9) with a
  # standard error of 1.1e4 on the log scale and said nothing
  # (dev/fams2-xbeta-noends.txt, the lane build before the check)
  set.seed(31)
  x <- rnorm(300)
  mu <- stats::plogis(0.2 + 0.5 * x)
  yb <- stats::rbeta(300, mu * 5, (1 - mu) * 5)
  d <- data.frame(x = x, y = yb)
  msg <- "this fit does not place it: kappa ran to 0 (at most"
  allow_warnings(frm(bf(y ~ x), family = xbeta(), data = d), msg,
                 require = msg)
  # PIN, the re-check's n3: kappa ~ x ran to 1.3e-10 at some rows and
  # 0.024 at others, with standard errors of 22 and 9.4, and the check,
  # which read the largest row, said nothing (dev/fams2-rev2-guards.txt)
  set.seed(31)
  x4 <- rnorm(400)
  mu4 <- stats::plogis(0.2 + 0.5 * x4)
  d4 <- data.frame(x = x4, y = stats::rbeta(400, mu4 * 5, (1 - mu4) * 5))
  msg <- "the standard error of kappa_(Intercept) is"
  allow_warnings(frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d4),
                 msg, require = msg)
  # PIN, the final check's K1: a real slope in log kappa runs some rows
  # to 2e-10 while the fit places kappa (standard errors 0.26 and 1.38,
  # 28 to 31 log-likelihood units above Beta()), and a rule on the
  # smallest row warned there (dev/fams2-rev3-kappa.txt, 3b)
  for (s in 2:3) {
    set.seed(300 + s)
    x5 <- stats::runif(2000, -4.7, 0)
    set.seed(400 + s)
    k5 <- exp(-1 + 3 * x5)
    z5 <- stats::rbeta(2000, 100, 100)
    d5 <- data.frame(x = x5, y = (1 + 2 * k5) * z5 - k5)
    expect_no_warning(frm(bf(y ~ 1, kappa ~ x), family = xbeta(), data = d5))
  }
  expect_no_warning(frm(bf(y ~ x, kappa = 0.1), family = xbeta(),
                        data = d))
  # one end is enough to place kappa
  d$y[1:10] <- 0
  expect_no_warning(frm(bf(y ~ x), family = xbeta(), data = d))
  # No end, and yet kappa is placed by the interior shape and grows: the
  # reviewer's A1 data (dev/fams2-rev-a1.R) reach kappa near 47 at phi
  # 2e5, 6 log-likelihood units above Beta(). The round-1 check warned
  # "kappa goes to 0 ... Use Beta()" there, which was false.
  skip_on_cran()
  set.seed(4101)
  z <- stats::rbeta(2000, 100, 100)
  da <- data.frame(y = 3 * z - 1)
  seen <- character()
  fa <- withCallingHandlers(frm(y ~ 1, family = xbeta(), data = da),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("kappa ran to 0", seen, fixed = TRUE)))
  fb <- frm(y ~ 1, family = Beta(), data = da)
  expect_gt(as.numeric(stats::logLik(fa)), as.numeric(stats::logLik(fb)))
  expect_gt(frmtmb:::eval_dpars(fa)[[1]]$kappa[1], 1)
  # PIN, the re-check's n3: at n 2000, seed 7, kappa stopped at 3.4e-5
  # on its way to 0, above the 1e-6 the fitted value is read against,
  # with a standard error of 96 on the log scale, and nothing was said
  # (dev/fams2-rev2-kappa9.txt). The standard error says it
  sim_beta <- function(n, phi, seed) {
    set.seed(seed)
    x <- rnorm(n)
    mu <- stats::plogis(0.2 + 0.5 * x)
    data.frame(x = x, y = stats::rbeta(n, mu * phi, (1 - mu) * phi))
  }
  msg <- "the standard error of kappa_(Intercept) is"
  allow_warnings(frm(y ~ x, family = xbeta(), data = sim_beta(2000, 5, 7)),
                 msg, require = msg)
  # and where the interior shape does place kappa (0.248, standard error
  # 1.35, 1.34 log-likelihood units above Beta()), nothing is said
  expect_no_warning(frm(y ~ x, family = xbeta(),
                        data = sim_beta(400, 50, 9)))
})

test_that("xbeta fits a precise response, with its ends", {
  skip_on_cran()
  # PIN: at phi 1e4 and above RTMB::dbeta()'s gradient, read with its
  # first argument on the tape, is NaN, and the round-1 build stopped at
  # "NA/NaN gradient evaluation" on this design (the reviewer's seed
  # 4200, dev/fams2-rev-phi5k.txt; dev/fams2-p1-pins.R replays it)
  ref <- function(y, mu, phi, kappa) sum(ref_xbeta(y, mu, phi, kappa))
  for (phi in c(1e4, 2e4)) {
    set.seed(4200)
    z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
    d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
    expect_gt(sum(d$y == 0), 0)
    # nlminb may call this false convergence at a precision in the
    # thousands; the gradient check below is what says it is the optimum
    fit <- allow_warnings(frm(y ~ 1, family = xbeta(), data = d),
                          "false convergence")
    dp <- frmtmb:::eval_dpars(fit)[[1]]
    ll <- as.numeric(stats::logLik(fit))
    r <- ref(d$y, dp$mu[1], dp$phi[1], dp$kappa[1])
    expect_lt(abs(ll - r), 1e-8 * abs(r), label = paste("phi", phi))
    expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
    # the estimate is the reference's own optimum: its gradient there,
    # by central differences on the link scale, is small against the
    # curvature the fit reports
    est <- c(stats::qlogis(dp$mu[1]), log(dp$phi[1]), log(dp$kappa[1]))
    g <- vapply(1:3, function(j) {
      h <- replace(numeric(3), j, 1e-5)
      f <- function(p) ref(d$y, stats::plogis(p[1]), exp(p[2]), exp(p[3]))
      (f(est + h) - f(est - h)) / 2e-5
    }, 0)
    se <- sqrt(diag(vcov(fit, full = TRUE)))[c("(Intercept)",
                                               "phi_(Intercept)",
                                               "kappa_(Intercept)")]
    # a unit of gradient times a standard error is a move of that many
    # standard errors in log-likelihood: well under one
    expect_lt(max(abs(g) * se), 1e-2, label = paste("phi", phi))
  }
  # kappa held at a constant sends the boundary rows through the
  # incomplete beta near its switch at phi 2e5, where the round-1 build
  # was off by 1.68 in logLik (dev/fams2-rev-kfix.R)
  set.seed(4200)
  z <- stats::rbeta(1000, 0.05 / 1.1 * 2e5, (1 - 0.05 / 1.1) * 2e5)
  d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
  fit <- frm(bf(y ~ 1, kappa = 0.05), family = xbeta(), data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  ll <- as.numeric(stats::logLik(fit))
  r <- ref(d$y, dp$mu[1], dp$phi[1], 0.05)
  expect_lt(abs(ll - r), 1e-8 * abs(r))
})

test_that("the incomplete beta has third derivatives where pbeta has none", {
  # RTMB::pbeta()'s third derivatives are NaN at this ordinary point,
  # which is why an xbeta fit with a random effect used to stop at
  # "NA/NaN gradient evaluation" (dev/fams2-pbeta-third.R)
  ours <- function(p) frmtmb:::log_ibeta_half(p[1], exp(p[2]), exp(p[3]))
  theirs <- function(p) log(RTMB::pbeta(p[1], exp(p[2]), exp(p[3])))
  p0 <- c(0.268941421, log(0.22313016), log(2.718282))
  mk <- function(f) {
    F <- RTMB::MakeTape(f, c(0.2, 0, 0))
    J <- F$jacfun()
    list(F = F, J = J, T3 = J$jacfun()$jacfun())
  }
  O <- mk(ours)
  expect_true(all(is.finite(O$T3(p0))))
  # FLIP this expectation, do not delete it, once RTMB's pbeta() has
  # finite third derivatives here: log_ibeta_half() could then use it
  expect_false(all(is.finite(mk(theirs)$T3(p0))))
  # x runs from the tail through m, the fraction's switch point, where
  # round 1 lost up to 1.0 in the log value at shape 1e6, and past it;
  # shapes to 8e5 (dev/reviews/2026-09-29-fams2.md, B2). Before punch 2
  # the error grew with the shapes, from the prefactor
  # a log x + b log(1 - x) - log B(a, b), whose terms cancel; formed
  # around the mean it holds 3e-13 to shapes of 1e7
  # (dev/fams2-p2-rev2-ibeta.txt), and the tolerance is 1.4e-12
  for (a in c(0.3, 4, 60, 3e3, 2e5)) for (b in c(0.5, 7, 300, 7e3, 8e5)) {
    s <- a + b
    m <- (a + 1) / (s + 2)
    sd <- sqrt(a * b / (s * s * (s + 1)))
    xs <- c(m, a / s, m * (1 - 1e-3), m * (1 + 1e-3), m - 2 * sd,
            m + 2 * sd, m - 5.5 * sd, m + 5.5 * sd, m - 6 * sd, m / 3, 0.45)
    for (x in xs[xs > 0 & xs < 0.5]) {
      r <- stats::pbeta(x, a, b, log.p = TRUE)
      if (!is.finite(r) || r < -600) next
      p <- c(x, log(a), log(b))
      expect_lt(abs(O$F(p) - r), 100 * ULPS * max(1, abs(r)),
                label = sprintf("a %g b %g x %.6g", a, b, x))
      expect_true(all(is.finite(O$T3(p))),
                  label = sprintf("third at a %g b %g x %.6g", a, b, x))
      # the x derivative is exactly dbeta / I (where it does not underflow),
      # relative and floored at one as the value is: 5.5 sd above m, where
      # I is near 1, it falls to 3.5e-3 and below and holds there to
      # 1.3e-12 absolutely, not relatively
      gx <- exp(stats::dbeta(x, a, b, log = TRUE) - r)
      if (!(gx > 0)) next
      expect_lt(abs(O$J(p)[1] - gx), 1e3 * ULPS * max(gx, 1),
                label = sprintf("d/dx at a %g b %g x %.6g", a, b, x))
    }
  }
  # PIN, the re-check's n1: RTMB::pbeta()'s derivatives are all NaN at
  # x == a / (a + b) exactly, which the pbeta branch read, so the
  # gradient here was NaN (dev/fams2-rev2-mean.txt)
  p <- c(0.3, log(300), log(700))
  expect_true(all(is.finite(O$J(p))))
  expect_true(all(is.finite(O$T3(p))))
  expect_lt(abs(O$F(p) - stats::pbeta(0.3, 300, 700, log.p = TRUE)),
            100 * ULPS)
  # and so an xbeta pair of rows (0, 0.3) at mu == q == 1/6 exactly
  fam <- frmtmb:::as_frmtmb_family(xbeta())
  X <- RTMB::MakeTape(function(p) {
    sum(fam$lpdf(c(0, 0.3), list(mu = p[1], phi = p[2], kappa = p[3]),
                 list()))
  }, c(0.2, 100, 0.1))
  expect_true(all(is.finite(X$jacfun()(c(1 / 6, 3000, 0.25)))))
  expect_true(all(is.finite(
    X$jacfun()$jacfun()$jacfun()(c(1 / 6, 3000, 0.25)))))
  # value and gradient join across x = m, where round 1's gradient read
  # 613 and 809 on the two sides at a = 3e4, b = 7e4 against 549
  a <- 3e4
  b <- 7e4
  m <- (a + 1) / (a + b + 2)
  gx <- exp(stats::dbeta(m, a, b, log = TRUE) -
              stats::pbeta(m, a, b, log.p = TRUE))
  for (x in c(m * (1 - 1e-12), m, m * (1 + 1e-12))) {
    p <- c(x, log(a), log(b))
    expect_lt(abs(O$J(p)[1] - gx), 1e-8 * gx)
    expect_true(all(is.finite(O$T3(p))))
  }
})

test_that("zero_inflated_beta_binomial's density is brms's", {
  fam <- zero_inflated_beta_binomial()
  y <- c(0, 0, 1, 5, 10)
  size <- c(10, 3, 10, 10, 10)
  dp <- list(mu = 0.35, phi = 4, zi = 0.25)
  ref <- ref_zibb(y, size, 0.35, 4, 0.25)
  got <- fam$lpdf(y, lapply(dp, rep, length.out = 5L), list(trials = size))
  expect_lt(rel_diff(got, ref), ULPS)
  tape <- c(dp, list(.eta_mu = stats::qlogis(0.35), .eta_phi = log(4),
                     .eta_zi = stats::qlogis(0.25)))
  got <- fam$lpdf(y, lapply(tape, rep, length.out = 5L), list(trials = size))
  expect_lt(rel_diff(got, ref), ULPS)
})

test_that("hurdle_cumulative's density is brms's under every link", {
  y <- c(0, 1, 2, 3, 4, 4)
  for (lk in names(hc_cdf)) {
    fam <- hurdle_cumulative(lk)
    tau <- c(-1, 0.3, 1.5)
    raw <- c(tau[1], log(diff(tau)))
    for (eta in c(-3, 0, 2)) {
      dp <- list(mu = rep(eta, 6), hu = rep(0.3, 6), disc = rep(1.4, 6))
      got <- fam$lpdf(y, dp, list(), list(tau_raw = raw))
      expect_lt(rel_diff(got, ref_hc(y, eta, 0.3, 1.4, tau, lk)), ULPS,
                label = paste(lk, eta))
    }
  }
})

# ------------------------------------------------------------ fit agreement

test_that("zero_inflated_beta_binomial matches glmmTMB's betabinomial", {
  skip_if_not_installed("glmmTMB")
  d <- sim_zibb(101, n = 600)
  fit <- frm(bf(y | trials(tr) ~ x, zi ~ x),
             family = zero_inflated_beta_binomial(), data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  ll <- as.numeric(stats::logLik(fit))
  expect_lt(abs(ll - sum(ref_zibb(d$y, d$tr, dp$mu, dp$phi, dp$zi))),
            ULPS * abs(ll))
  ref <- glmmTMB::glmmTMB(cbind(y, tr - y) ~ x, ziformula = ~ x,
                          family = glmmTMB::betabinomial(), data = d)
  expect_lt(abs(ll - as.numeric(stats::logLik(ref))), 1e-8 * abs(ll))
  fe <- fixef_by_dpar(fit)
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  expect_lt(max(abs(fe$mu - glmmTMB::fixef(ref)$cond) /
                  se[c("(Intercept)", "x")]), 1e-2)
  expect_lt(max(abs(fe$zi - glmmTMB::fixef(ref)$zi) /
                  se[c("zi_(Intercept)", "zi_x")]), 1e-2)
})

test_that("hurdle_cumulative factorizes into a glm and a polr", {
  # IDENTITY: with separate predictors the log-likelihood is a bernoulli
  # on 1{y = 0} and a cumulative model over the rows above zero, sharing
  # no parameter, so the joint ML fit is the two separate fits
  d <- sim_hc(102, n = 600)
  fit <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(), data = d)
  r_hu <- stats::glm(I(y == 0) ~ z, family = stats::binomial, data = d)
  r_ord <- MASS::polr(factor(y, ordered = TRUE) ~ x, data = d[d$y > 0, ])
  ll <- as.numeric(stats::logLik(fit))
  ll2 <- as.numeric(stats::logLik(r_hu)) + as.numeric(stats::logLik(r_ord))
  expect_lt(abs(ll - ll2), 1e-8 * abs(ll))
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  tau <- frmtmb:::ord_tau_from_raw(fit$estimates$tau_raw, TRUE)
  expect_lt(abs(ll - sum(ref_hc(d$y, dp$mu, dp$hu, 1, tau, "logit"))),
            ULPS * abs(ll))
  fx <- fixef(fit)
  expect_lt(max(abs(fx[c("Intercept[1]", "Intercept[2]", "Intercept[3]"),
                       "Estimate"] - r_ord$zeta) /
                  fx[c("Intercept[1]", "Intercept[2]", "Intercept[3]"),
                     "Est.Error"]), 1e-2)
  expect_lt(abs(fx["hu_z", "Estimate"] - stats::coef(r_hu)[["z"]]) /
              fx["hu_z", "Est.Error"], 1e-2)
})

test_that("xbeta and the zi beta-binomial fit with a random intercept", {
  skip_on_cran()
  # PIN: xbeta on RTMB::pbeta() stopped at "NA/NaN gradient evaluation"
  # on this design and on each of ten seeds (dev/fams2-validate.R, 2)
  d <- sim_xbeta(103, n = 600, ngrp = 20)
  expect_gt(sum(d$y == 0) + sum(d$y == 1), 20)
  fit <- frm(bf(y ~ x + (1 | g)), family = xbeta(), data = d)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  skip_if_not_installed("glmmTMB")
  dz <- sim_zibb(104, n = 600, ngrp = 20)
  fz <- frm(bf(y | trials(tr) ~ x + (1 | g), zi ~ x),
            family = zero_inflated_beta_binomial(), data = dz)
  ref <- glmmTMB::glmmTMB(cbind(y, tr - y) ~ x + (1 | g), ziformula = ~ x,
                          family = glmmTMB::betabinomial(), data = dz)
  ll <- as.numeric(stats::logLik(fz))
  expect_lt(abs(ll - as.numeric(stats::logLik(ref))), 1e-7 * abs(ll))
})

# ---------------------------------------------------------------- post-fit

test_that("fitted() is brms's posterior_epred for all three families", {
  d <- sim_xbeta(5, n = 200)
  fit <- frm(bf(y ~ x), family = xbeta(), data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  # brms:::posterior_epred_xbeta(), written out
  a <- dp$mu * dp$phi
  b <- (1 - dp$mu) * dp$phi
  dd <- 1 + 2 * dp$kappa
  q0 <- dp$kappa / dd
  q1 <- (1 + dp$kappa) / dd
  t3 <- stats::pbeta(q1, a, b)
  m <- 1 + dd * dp$mu * (stats::pbeta(q1, a + 1, b) -
                           stats::pbeta(q0, a + 1, b)) -
    dp$kappa * (t3 - stats::pbeta(q0, a, b)) - t3
  expect_lt(rel_diff(fitted(fit)[, "Estimate"], m), ULPS)
  expect_lt(rel_diff(fitted(fit, dpar = "kappa")[, "Estimate"],
                     rep(dp$kappa, length.out = 200)), ULPS)

  dz <- sim_zibb(6, n = 200)
  fz <- frm(bf(y | trials(tr) ~ x), family = zero_inflated_beta_binomial(),
            data = dz)
  dp <- frmtmb:::eval_dpars(fz)[[1]]
  expect_lt(rel_diff(fitted(fz)[, "Estimate"],
                     dp$mu * dz$tr * (1 - dp$zi)), ULPS)

  dh <- sim_hc(7, n = 200)
  fh <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(), data = dh)
  dp <- frmtmb:::eval_dpars(fh)[[1]]
  tau <- frmtmb:::ord_tau_from_raw(fh$estimates$tau_raw, TRUE)
  # brms:::posterior_epred_hurdle_cumulative(): cbind(hu, (1 - hu) *
  # dcumulative())
  Fm <- stats::plogis(outer(-dp$mu, tau, "+"))
  e <- cbind(dp$hu, (cbind(Fm, 1) - cbind(0, Fm)) * (1 - dp$hu))
  P <- frm_linpred(fh, type = "response")
  expect_identical(colnames(P), as.character(0:4))
  expect_lt(max(abs(unname(P) - e) / e), ULPS)
  expect_lt(rel_diff(as.numeric(frm_linpred(fh, type = "zprob")), dp$hu),
            ULPS)
})

test_that("the variance functions are the densities' own second moments", {
  fx <- xbeta()
  dp <- list(mu = 0.4, phi = 5, kappa = 0.2)
  atoms <- exp(ref_xbeta(c(0, 1), 0.4, 5, 0.2))
  # each moment with the quadrature's own error estimate, which is what
  # the comparisons below are held to
  mom <- function(k) {
    q <- stats::integrate(function(y) y^k * exp(ref_xbeta(y, 0.4, 5, 0.2)),
                          0, 1, rel.tol = 1e-12)
    c(value = q$value + atoms[[2]], err = q$abs.error)
  }
  m0 <- mom(0)
  m1 <- mom(1)
  m2 <- mom(2)
  # the density and its two atoms carry all the mass
  expect_lt(abs(m0[["value"]] + atoms[1] - 1), 10 * m0[["err"]] + ULPS)
  expect_lt(abs(fx$post$mean_fn(dp, list()) - m1[["value"]]),
            10 * m1[["err"]] + ULPS * m1[["value"]])
  v <- m2[["value"]] - m1[["value"]]^2
  expect_lt(abs(fx$post$var_fn(dp, list()) - v),
            10 * (m2[["err"]] + 2 * m1[["err"]]) + ULPS * v)
  fz <- zero_inflated_beta_binomial()
  dp <- list(mu = 0.3, phi = 3, zi = 0.2)
  y <- 0:12
  p <- exp(ref_zibb(y, 12, 0.3, 3, 0.2))
  expect_lt(abs(sum(p) - 1), ULPS * length(p))
  m1 <- sum(y * p)
  expect_lt(rel_diff(fz$post$mean_fn(dp, list(trials = 12)), m1),
            ULPS * length(p))
  expect_lt(rel_diff(fz$post$var_fn(dp, list(trials = 12)),
                     sum(y^2 * p) - m1^2), 1e3 * ULPS)
})
test_that("beta_binomial and its zero-inflated form give pearson residuals", {
  # PIN: beta_binomial() had no variance function, so pearson residuals
  # were refused with "has no variance function"; seen failing on the
  # base build (dev/fams2-findings.md)
  d <- sim_zibb(8, n = 150)
  fb <- frm(bf(y | trials(tr) ~ x), family = beta_binomial(), data = d)
  r <- residuals(fb, type = "pearson")
  dp <- frmtmb:::eval_dpars(fb)[[1]]
  v <- d$tr * dp$mu * (1 - dp$mu) * (dp$phi + d$tr) / (dp$phi + 1)
  expect_lt(rel_diff(r[, "Estimate"], (d$y - d$tr * dp$mu) / sqrt(v)),
            1e3 * ULPS)
  fz <- frm(bf(y | trials(tr) ~ x), family = zero_inflated_beta_binomial(),
            data = d)
  expect_identical(nrow(residuals(fz, type = "pearson")), 150L)
})

test_that("a hurdle ordinal fit reads codes 0..K through every method", {
  d <- sim_hc(9, n = 300)
  lv <- c("none", "low", "mid", "high", "top")
  d$f <- factor(lv[d$y + 1L], levels = lv, ordered = TRUE)
  fi <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(), data = d)
  ff <- frm(bf(f ~ x, hu ~ z), family = hurdle_cumulative(), data = d)
  # the first level is the hurdle, so the two fits are one model
  expect_lt(abs(as.numeric(stats::logLik(ff)) -
                  as.numeric(stats::logLik(fi))),
            ULPS * abs(as.numeric(stats::logLik(fi))))
  expect_identical(colnames(frm_linpred(ff, type = "response")), lv)
  s <- simulate(ff, nsim = 1, seed = 1)[[1]]
  expect_true(is.ordered(s))
  expect_identical(levels(s), lv)
  expect_true(all(simulate(fi, nsim = 1, seed = 1)[[1]] %in% 0:4))
  set.seed(2)
  pr <- predict(fi, ndraws = 200)
  expect_identical(dim(pr), c(300L, 5L))
  # the expected category is scored by the codes, so the hurdle scores 0
  nd <- data.frame(x = 0.5, z = -0.2)
  P <- frm_linpred(fi, newdata = nd, type = "response")
  ce <- allow_warnings(
    conditional_effects(fi, effects = "x", categorical = FALSE,
                        int_conditions = list(x = 0.5),
                        conditions = data.frame(z = -0.2)),
    "Predictions are treated as continuous")
  expect_lt(abs(ce[[1]]$estimate__[1] - sum(P * 0:4)), ULPS * sum(P * 0:4))
  # P(Y = 0) is hu alone, so its band is hu's own delta-method band
  cc <- conditional_effects(fi, effects = "x", categorical = TRUE,
                            int_conditions = list(x = 0.5),
                            conditions = data.frame(z = -0.2))[[1]]
  hu <- fitted(fi, newdata = nd, dpar = "hu")
  row0 <- cc[cc$cats__ == "0", ]
  expect_lt(abs(row0$estimate__ - hu[, "Estimate"]), ULPS * hu[, "Estimate"])
  expect_lt(abs(row0$se__ - hu[, "Est.Error"]) / hu[, "Est.Error"], 1e-4)
})

test_that("thres(x = ) counts the categories above the hurdle", {
  d <- sim_hc(10, n = 300)
  fit <- allow_warnings(
    frm(bf(y | thres(5) ~ x, hu ~ z), family = hurdle_cumulative(),
        data = d),
    "response categories that no row takes",
    require = "response categories that no row takes")
  expect_identical(ncol(frm_linpred(fit, type = "response")), 7L)
})

test_that("disc is held at one unless the formula models it", {
  d <- sim_hc(11, n = 300)
  f1 <- frm(bf(y ~ x), family = hurdle_cumulative(), data = d)
  expect_identical(summary(f1)$fixed_dpars[["disc"]], 1)
  f2 <- frm(bf(y ~ x, disc ~ 0 + z), family = hurdle_cumulative("probit"),
            data = d)
  dp <- frmtmb:::eval_dpars(f2)[[1]]
  tau <- frmtmb:::ord_tau_from_raw(f2$estimates$tau_raw, TRUE)
  ll <- as.numeric(stats::logLik(f2))
  expect_lt(abs(ll - sum(ref_hc(d$y, dp$mu, dp$hu, dp$disc, tau,
                                "probit"))), ULPS * abs(ll))
  expect_match(frmtmb:::family_link_str(family(f2)),
               "cdf = probit; hu = logit; disc = log", fixed = TRUE)
})

test_that("an intercept in disc warns unless a prior holds it", {
  # PIN: disc ~ 1 + z fitted silently with standard errors of 12 to 67
  # on the thresholds and coefficients (dev/fams2-rev-disc1.txt). brms
  # accepts it and identifies it by a normal(0, 1) default prior
  d <- sim_hc(12, n = 400)
  msg <- "have no usable standard error"
  allow_warnings(frm(bf(y ~ x, disc ~ 1 + z), family = hurdle_cumulative(),
                     data = d), msg, require = msg)
  # PIN, the re-check's n4: a prior on the thresholds pins the scale
  # through that prior, yet the warning said "no usable standard error"
  # where the largest was 1.71 (dev/fams2-rev2-guards.txt, 4)
  msg <- "only the priors on the thresholds or on the coefficients of mu"
  allow_warnings(frm(bf(y ~ x, disc ~ 1 + z), family = hurdle_cumulative(),
                     data = d,
                     prior = set_prior("normal(0, 3)", class = "Intercept")),
                 msg, require = msg)
  expect_no_warning(frm(bf(y ~ x, disc ~ 0 + z),
                        family = hurdle_cumulative(), data = d))
  expect_no_warning(fp <- frm(
    bf(y ~ x, disc ~ 1 + z), family = hurdle_cumulative(), data = d,
    prior = set_prior("normal(0, 1)", class = "Intercept", dpar = "disc")))
  expect_true(all(is.finite(sqrt(diag(vcov(fp))))))
})

test_that("conditional_effects() and emmeans() run on the three families", {
  skip_if_not_installed("emmeans")
  fx <- frm(bf(y ~ x), family = xbeta(), data = sim_xbeta(12, n = 200))
  ce <- conditional_effects(fx)
  expect_true(all(ce[[1]]$estimate__ > 0 & ce[[1]]$estimate__ < 1))
  expect_s4_class(emmeans::emmeans(fx, ~ x), "emmGrid")
  fz <- frm(bf(y | trials(tr) ~ x), family = zero_inflated_beta_binomial(),
            data = sim_zibb(13, n = 200))
  ce <- conditional_effects(fz, conditions = data.frame(tr = 10))
  expect_true(all(ce[[1]]$estimate__ > 0 & ce[[1]]$estimate__ < 10))
  fh <- frm(bf(y ~ x, hu ~ z), family = hurdle_cumulative(),
            data = sim_hc(14, n = 200))
  expect_s4_class(emmeans::emmeans(fh, ~ x), "emmGrid")
  ce <- conditional_effects(fh, dpar = "hu", effects = "z")
  expect_true(all(ce[[1]]$estimate__ > 0 & ce[[1]]$estimate__ < 1))
})

test_that("cse() is cs() under brms's second name", {
  # PIN: the base build had no cse() and failed on it
  skip_if_not_installed("brms")
  inhaler <- brms::inhaler
  f1 <- frm(rating ~ treat + period + cs(carry), data = inhaler,
            family = sratio())
  f2 <- frm(rating ~ treat + period + cse(carry), data = inhaler,
            family = sratio())
  expect_identical(as.numeric(stats::logLik(f2)),
                   as.numeric(stats::logLik(f1)))
  expect_identical(fixef(f2), fixef(f1))
  expect_error(frm(rating ~ treat * cse(carry), data = inhaler,
                   family = sratio()), "cse() makes a whole term",
               fixed = TRUE)
})

test_that("the compatibility registry lists the three families", {
  ft <- frm_compat_features()
  fams <- c("xbeta", "zero_inflated_beta_binomial", "hurdle_cumulative")
  expect_true(all(fams %in% ft$key[ft$kind == "family"]))
  for (f in fams) {
    expect_identical(frm_compat(f, "simulate")$status, "works", info = f)
    expect_identical(frm_compat(f, "residuals_osa")$status, "refused",
                     info = f)
  }
  expect_identical(frm_compat("zero_inflated_beta_binomial",
                              "trials()")$status, "works")
  expect_identical(frm_compat("hurdle_cumulative", "fitted")$status,
                   "conditional")
})
