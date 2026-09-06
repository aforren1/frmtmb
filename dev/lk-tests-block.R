
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

test_that("the links with no exact robust form leave the field absent", {
  # softplus: log(log1p(exp(eta))) is exact to eta = -745, where the
  # softplus underflows, and no branch-free form recovers it
  sp <- frmtmb:::frmtmb_links$softplus
  expect_null(sp$log_eta)
  expect_equal(log(sp$linkinv(-700)), -700, tolerance = 1e-12)
  expect_identical(log(sp$linkinv(-746)), -Inf)
  # sqrt and `1/mu^2` never lose a digit through the round trip
  expect_null(frmtmb:::frmtmb_links$sqrt$log_eta)
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
      expect_robust_at(bernoulli(lk), y01, list(), list(mu = 0), "mu", e0)
      expect_robust_at(frmtmb:::fam_binomial(lk), ybin, list(trials = 5),
                       list(mu = 0), "mu", e0)
      expect_robust_at(Beta(lk), c(0.2, 0.5, 0.9), list(),
                       list(mu = 0, phi = log(5)), "mu", e0)
    }
  }
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
  d$yp <- stats::rpois(n, (1.2 + 0.4 * d$x)^2)
  fp <- suppressWarnings(frm(bf(yp ~ x),
                             family = stats::poisson(link = "sqrt"),
                             data = d, control = ctl))
  gp <- stats::glm(yp ~ x, data = d,
                   family = stats::poisson(link = "sqrt"), control = gctl)
  expect_lt(max(abs(fixef(fp)$mu - stats::coef(gp))), 1e-6)

  # `1/mu^2` is the inverse Gaussian canonical link and stats spells it
  # the same way brms does, so this one comparison is available
  d$yi <- stats::rgamma(n, 8, 8 / exp(0.4 + 0.2 * d$x))
  fi <- suppressWarnings(frm(bf(yi ~ x),
                             family = stats::inverse.gaussian(),
                             data = d, control = ctl))
  expect_identical(family(fi)$links$mu$name, "1/mu^2")
  gi <- suppressWarnings(stats::glm(
    yi ~ x, data = d, family = stats::inverse.gaussian(),
    start = c(1 / mean(d$yi)^2, 0), control = gctl))
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
