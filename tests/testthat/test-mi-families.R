# mi() on a response of any family brms 2.23.0 takes it for. A value no
# density reads leaves frm()'s objective, so the fit is the complete-case one,
# exactly; a value a density reads stays latent, on a chart that keeps it in
# the family's support. The identities against brms's log_prob are in
# test-mi-mm-brms-l1.R (gated).

mi_data <- function(fam, seed = 169, n = 150, nmiss = 12) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  eta <- 0.5 + 0.4 * d$x
  d$y <- switch(fam,
    gaussian = rnorm(n, eta, 1),
    student = eta + rt(n, 4),
    skew_normal = eta + rnorm(n) + abs(rnorm(n)),
    beta = rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8),
    Gamma = rgamma(n, 3, 3 / exp(eta)),
    weibull = rweibull(n, 2, exp(eta)),
    exponential = rexp(n, 1 / exp(eta)),
    frechet = 1 / rweibull(n, 3, exp(-eta)),
    inverse.gaussian = rgamma(n, 4, 4 / exp(eta)),
    lognormal = rlnorm(n, eta, 0.6),
    exgaussian = eta + rnorm(n, 0, 0.5) + rexp(n, 1),
    gen_extreme_value = eta + 0.6 * (rexp(n)^(-0.1) - 1) / 0.1,
    von_mises = ((2 * atan(0.5 * d$x) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) -
      pi,
    asym_laplace = eta + rexp(n) - rexp(n),
    poisson = rpois(n, exp(eta)))
  d$full <- d$y
  d$y[sample(n, nmiss)] <- NA
  d
}

mi_family <- function(fam) {
  switch(fam,
    beta = Beta(),
    Gamma = stats::Gamma(link = "log"),
    inverse.gaussian = stats::inverse.gaussian(link = "log"),
    get(fam, asNamespace("frmtmb"))())
}

mi_fams <- c("gaussian", "student", "skew_normal", "beta", "Gamma", "weibull",
             "exponential", "frechet", "inverse.gaussian", "lognormal",
             "exgaussian", "gen_extreme_value", "von_mises")

# the fit of a family with a kink warns about its gradient on both arms
mi_kink_warn <- c("Optimizer did not report convergence",
                  "Large maximum absolute gradient")

test_that("a value no density reads is the complete-case fit, exactly", {
  for (fam in c(mi_fams, "asym_laplace")) {
    dd <- mi_data(fam)
    a <- allow_warnings(frm(bf(y | mi() ~ x), family = mi_family(fam),
                            data = dd), mi_kink_warn)
    b <- allow_warnings(frm(bf(y ~ x), family = mi_family(fam),
                            data = dd[!is.na(dd$y), ]), mi_kink_warn)
    la <- as.numeric(logLik(a))
    lb <- as.numeric(logLik(b))
    expect_lt(abs(la - lb) / abs(lb), 1e-7, label = fam)
    expect_lt(abs(AIC(a) - AIC(b)) / abs(AIC(b)), 1e-7, label = fam)
    expect_lt(abs(BIC(a) - BIC(b)) / abs(BIC(b)), 1e-7, label = fam)
    # the rows stay in the frame, and predict like any other
    expect_null(a$estimates$miss)
    ft <- fitted(a)
    expect_identical(nrow(ft), nrow(dd))
    expect_true(all(is.finite(ft[, "Estimate"])), info = fam)
    # the coefficient vector, in its own size: a kink has no standard error
    ba <- frmtmb:::fixef_by_dpar(a)$mu
    bb <- frmtmb:::fixef_by_dpar(b)$mu
    expect_lt(max(abs(ba - bb)) / max(abs(bb)), 1e-3, label = fam)
  }
})

test_that("a student fit reads no Laplace error any more", {
  # 0.78.0 fitted these rows latent, -0.165 per missing value at nu near 4,
  # and the fit moved toward a gaussian as they grew
  dd <- mi_data("student", nmiss = 60)
  a <- frm(bf(y | mi() ~ x), family = student(), data = dd)
  b <- frm(bf(y ~ x), family = student(), data = dd[!is.na(dd$y), ])
  ea <- unlist(a$estimates[c("beta", "betad")])
  eb <- unlist(b$estimates[c("beta", "betad")])
  expect_lt(max(abs(ea - eb)) / max(abs(eb)), 1e-5)
})

test_that("a discrete or inflated family is refused in brms's words", {
  dd <- mi_data("poisson")
  expect_error(frm(bf(y | mi() ~ x), family = poisson(), data = dd,
                   dry_run = "frame"),
               "Argument 'mi' is not supported for family 'poisson'",
               fixed = TRUE)
  for (fam in c("hurdle_gamma", "zero_inflated_beta", "shifted_lognormal",
                "xbeta")) {
    d2 <- mi_data("beta")
    d2$y[!is.na(d2$y)] <- pmin(d2$y[!is.na(d2$y)], 0.9)
    expect_error(frm(bf(y | mi() ~ x), family = get(fam)(), data = d2,
                     dry_run = "frame"),
                 "Argument 'mi' is not supported", fixed = TRUE,
                 info = fam)
  }
})

test_that("a response with no observed value is refused by name", {
  dd <- mi_data("Gamma")
  dd$y <- NA_real_
  for (fm in list(stats::Gamma(link = "log"), gaussian())) {
    expect_error(frm(bf(y | mi() ~ x), family = fm, data = dd,
                     dry_run = "frame"),
                 "has no observed value", fixed = TRUE)
  }
})

mi_read_data <- function(fam, seed = 176, nmiss = 12) {
  dd <- mi_data(fam, seed = seed, nmiss = nmiss)
  dd$resp <- dd$z + 0.5 * dd$full + rnorm(nrow(dd), 0, 0.5)
  dd
}

mi_read_fit <- function(fam) {
  allow_warnings(
    frm(mvbf(bf(resp ~ mi(y) + z), bf(y | mi() ~ x, family = mi_family(fam)),
             rescor = FALSE), data = mi_read_data(fam)),
    "Laplace approximation of the")
}

test_that("a value a density reads stays latent, on the family's chart", {
  for (fam in mi_fams) {
    fit <- mi_read_fit(fam)
    expect_true(is.finite(as.numeric(logLik(fit))), info = fam)
    expect_length(fit$estimates$miss, 12L)
    expect_length(fit$frame$mi_map$y$unread, 0L)
  }
  fg <- mi_read_fit("Gamma")
  fb <- mi_read_fit("beta")
  fn <- mi_read_fit("gaussian")
  expect_identical(fg$frame$mi_map$y$chart, "log")
  expect_identical(fb$frame$mi_map$y$chart, "logit")
  expect_identical(fn$frame$mi_map$y$chart, "identity")
  v <- frmtmb:::mi_values(fg, "y")
  expect_true(all(v > 0))
  w <- frmtmb:::mi_values(fb, "y")
  expect_true(all(w > 0 & w < 1))
  expect_true(all(is.finite(log(w))))
  # the observed rows are the data, unmoved
  dd <- mi_read_data("Gamma")
  expect_identical(v[!is.na(dd$y)], dd$y[!is.na(dd$y)])
})

test_that("only the rows a density reads stay latent", {
  # mi(y, idx = ) reads the rows of y its index names, so a response read in
  # part keeps those rows and lets the others go
  dd <- mi_data("Gamma", nmiss = 0L)
  dd$id <- seq_len(nrow(dd))
  dd$y[c(3L, 7L, 11L, 20L)] <- NA
  dd$read <- c(rep(TRUE, 12), rep(FALSE, nrow(dd) - 12))
  fit <- frm(mvbf(bf(y | mi() + index(id) ~ x,
                     family = stats::Gamma(link = "log")),
                  bf(z | subset(read) ~ mi(y, idx = id)),
                  rescor = FALSE), data = dd)
  m <- fit$frame$mi_map$y
  expect_setequal(m$rows, c(3L, 7L, 11L))
  expect_identical(m$unread, 20L)
  expect_length(fit$estimates$miss, 3L)
})

test_that("the chart's Jacobian is on the taped objective", {
  fit <- mi_read_fit("exponential")
  env <- fit$obj$env
  m <- fit$frame$mi_map$y
  par <- env$last.par.best
  mi_i <- which(names(par) == "miss")
  # the tape against the density by hand, in one latent value: the value's
  # exponential density and its reader's normal, with the chart's u
  bx <- frmtmb:::fixef_by_dpar(fit)
  dd <- mi_read_data("exponential")
  i <- m$rows[1L]
  u0 <- par[mi_i[1L]]
  eta_y <- bx$y_mu[[1L]] + bx$y_mu[[2L]] * dd$x[i]
  sg <- exp(bx$resp_sigma[[1L]])
  mz <- bx$resp_mu[["(Intercept)"]] + bx$resp_mu[["z"]] * dd$z[i]
  hand <- function(v) {
    stats::dexp(exp(v), exp(-eta_y), log = TRUE) + v +
      stats::dnorm(dd$resp[i], mz + bx$resp_mu[["miy"]] * exp(v), sg,
                   log = TRUE)
  }
  p2 <- par
  p2[mi_i[1L]] <- u0 + 0.3
  ta <- -(env$f(p2) - env$f(par))
  hd <- hand(u0 + 0.3) - hand(u0)
  expect_lt(abs(ta - hd) / abs(hd), 1e-8)
})

test_that("Laplace on a read value is exact where the integrand is normal", {
  # gaussian under a gaussian reader; elsewhere the error is the
  # integrand's, as NEWS states with its numbers
  fit <- mi_read_fit("gaussian")
  env <- fit$obj$env
  par <- env$last.par.best
  mi_i <- which(names(par) == "miss")
  base <- env$f(par)
  u0 <- par[mi_i[1L]]
  w <- 1 / sqrt((env$f(replace(par, mi_i[1L], u0 + 1e-3)) +
                   env$f(replace(par, mi_i[1L], u0 - 1e-3)) - 2 * base) / 1e-6)
  h <- function(u) {
    vapply(u, function(ui) exp(base - env$f(replace(par, mi_i[1L], ui))), 0)
  }
  ex <- stats::integrate(h, u0 - 40 * w, u0 + 40 * w, rel.tol = 1e-9,
                         subdivisions = 2000L)$value
  expect_lt(abs(sqrt(2 * pi) * w / ex - 1), 1e-5)
})

test_that("mi(sd) measures a charted latent on its own scale", {
  dd <- mi_data("lognormal")
  dd$y[is.na(dd$y)] <- dd$full[is.na(dd$y)]
  dd$s <- 0.1
  fit <- frm(bf(y | mi(s) ~ x), family = lognormal(), data = dd)
  expect_length(fit$estimates$miss, nrow(dd))
  expect_identical(fit$frame$mi_map$y$chart, "log")
  v <- frmtmb:::mi_values(fit, "y")
  expect_true(all(v > 0))
  # every latent sits near its measurement, at sd 0.1 on the data's scale
  expect_lt(stats::median(abs(v - dd$y)) / stats::sd(dd$y), 0.3)
})

test_that("a von_mises measurement at the end of the window fits", {
  set.seed(170)
  dd <- data.frame(x = rnorm(80))
  dd$y <- ((2 * atan(0.5 * dd$x) + rnorm(80, 0, 0.5) + pi) %% (2 * pi)) - pi
  dd$y[1L] <- pi
  dd$s <- 0.3
  fit <- allow_warnings(frm(bf(y | mi(s) ~ x), family = von_mises(), data = dd),
                        "Laplace approximation of the")
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("frm() and the sampler hold a von_mises value in its window", {
  set.seed(201)
  n <- 160
  dd <- data.frame(z = rnorm(n))
  ang <- rnorm(n, 0, 0.3)
  dd$y <- ((ang + pi) %% (2 * pi)) - pi
  dd$z <- dd$z + ang
  dd$y[1:12] <- NA
  f <- mvbf(bf(z ~ mi(y)), bf(y | mi() ~ 1, family = von_mises()),
            rescor = FALSE)
  fit <- allow_warnings(frm(f, data = dd), "Laplace approximation of the")
  expect_identical(fit$frame$mi_map$y$chart, "circle")
  expect_true(all(abs(frmtmb:::mi_values(fit, "y")[1:12]) < pi))
  obj <- frm(f, data = dd, dry_run = "frame")
  expect_identical(obj$mi_map$y$chart, "circle")
})

test_that("a charted mi() variable reaches its predictor on the data scale", {
  fit <- mi_read_fit("Gamma")
  dd <- mi_read_data("Gamma")
  b <- frmtmb:::fixef_by_dpar(fit)$resp_mu
  ft <- fitted(fit, resp = "resp")
  v <- frmtmb:::mi_values(fit, "y")
  eta <- b[["(Intercept)"]] + b[["miy"]] * v + b[["z"]] * dd$z
  expect_lt(max(abs(ft[, "Estimate"] - eta)) / max(abs(eta)), 1e-8)
})

test_that("a read value with no mode to expand about is refused", {
  dd <- mi_read_data("asym_laplace")
  f <- mvbf(bf(resp ~ mi(y) + z), bf(y | mi() ~ x, family = asym_laplace()),
            rescor = FALSE)
  expect_error(frm(f, data = dd), "no curvature at its kink", fixed = TRUE)
  # the objective a sampler reads is built
  expect_no_error(frm(f, data = dd, dry_run = "objective"))
})

test_that("the bounded charts stay finite at both ends of the tape", {
  # a value that rounds to the end of its interval is log(0) in the density,
  # which stopped beta at a mean near 1 and not near 0
  e <- frmtmb:::mi_chart_inv(c(-800, -40, 40, 800), "logit")
  expect_true(all(e > 0 & e < 1))
  expect_true(all(is.finite(stats::dbeta(e, 2, 0.5, log = TRUE))))
  expect_true(all(is.finite(frmtmb:::mi_chart_logjac(c(-800, 800), "logit"))))
  tp <- RTMB::MakeTape(function(p) {
    y <- frmtmb:::mi_chart_inv(p[1], "logit")
    log(y) - 0.5 * log1p(-y) + frmtmb:::mi_chart_logjac(p[1], "logit")
  }, 0)
  expect_true(is.finite(tp(40)))
  expect_true(is.finite(tp(-40)))
  # and it inverts where a measurement can be
  expect_true(all(is.finite(frmtmb:::mi_chart_fun(c(0, 1), "logit"))))
  expect_true(all(is.finite(frmtmb:::mi_chart_fun(c(-pi, pi), "circle"))))
})

test_that("a beta mean near 1 fits a read mi(), as a mean near 0 does", {
  # dev/refuse-f1-m4.R reads the whole grid of the review's 24 fits
  for (eta0 in c(3, -3)) {
    set.seed(200)
    n <- 200
    dd <- data.frame(x = rnorm(n))
    m <- stats::plogis(eta0 + 0.5 * dd$x)
    dd$y <- pmin(pmax(rbeta(n, m * 8, (1 - m) * 8), 1e-10), 1 - 1e-10)
    dd$z <- 0.5 + stats::qlogis(dd$y) / 5 + rnorm(n, 0, 0.5)
    dd$y[seq(5, 195, 10)] <- NA
    fit <- allow_warnings(
      frm(mvbf(bf(z ~ mi(y)), bf(y | mi() ~ x, family = Beta()),
               rescor = FALSE), data = dd),
      "Laplace approximation of the")
    expect_true(is.finite(as.numeric(logLik(fit))), info = eta0)
    expect_identical(fit$opt$convergence, 0L)
  }
})

test_that("a family's mi_chart is checked by name", {
  expect_error(
    frmtmb_family("x", dpars = "mu", links = list(mu = "identity"),
                  lpdf = function(y, dpars, aterms) y, mi_chart = "c"),
    "mi_chart", fixed = TRUE)
})

test_that("a family of one's own with bounded support fits unread rows", {
  # its valid_y reads the column, which holds no 0 at a missing row
  fam <- frmtmb_family(
    "tpos", dpars = c("mu", "sigma"),
    links = list(mu = "log", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dgamma(y, 1 / dpars[["sigma"]]^2,
                   dpars[["mu"]] * dpars[["sigma"]]^2, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0)) stop("tpos: response must be strictly positive")
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) 0.5),
    accepts_aterms = "mi")
  dd <- mi_data("Gamma")
  fit <- frm(bf(y | mi() ~ x), family = fam, data = dd)
  expect_true(is.finite(as.numeric(logLik(fit))))
})
