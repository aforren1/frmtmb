# A nonlinear fit from its default start reaches the maximum where the
# cause of a short stop is a sign pattern or a scale frm() can see:
# nl_sign_probe() and nl_par_units() (R/par-template.R). Each defect test
# below was seen to fail on 0.79.1 (dev/nlstart-findings.md).

sat_data <- function(seed) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z + rnorm(n, 0, 0.5)
  d
}
sat_body <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)

# every message and warning a call raises, by class, muffled
caught <- function(expr) {
  got <- list()
  val <- withCallingHandlers(expr,
    warning = function(w) {
      got[[length(got) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      got[[length(got) + 1L]] <<- m
      invokeRestart("muffleMessage")
    })
  list(value = val, conds = got)
}
has_class <- function(conds, cls) {
  any(vapply(conds, function(cnd) inherits(cnd, cls), NA))
}

test_that("a body whose maximum has both signs negative is reached", {
  for (seed in c(101, 201)) {
    d <- sat_data(seed)
    got <- caught(frm(sat_body, data = d))
    fit <- got$value
    # the maximum, from a start inside its basin
    ref <- caught(frm(sat_body, data = d, start = list(beta = c(
      "b1_(Intercept)" = -2, "b2_(Intercept)" = -0.5))))$value
    ll <- as.numeric(logLik(fit))
    ll_ref <- as.numeric(logLik(ref))
    # 0.79.1 stopped at the linear limit, about 70 below, at code 0
    expect_gte(ll, ll_ref - 1e-8 * abs(ll_ref))
    expect_lt(fit$estimates$beta[["b2_(Intercept)"]], 0)
    expect_true(has_class(got$conds, "frmtmb_nl_multistart"))
    # the stop it left warned on seed 201 (false convergence); that
    # warning, and the message of the restart it kept, belong to a fit
    # frm() does not return
    expect_false(any(vapply(got$conds, inherits, NA, "warning")))
    expect_false(has_class(got$conds, "frmtmb_nl_restart"))
  }
})

# the review's two saturating terms (dev/review-nlstart-fresh.R, body
# two_sat_rev, bernoulli, seed 33): every draw in the order it made them
two_sat_data <- function() {
  set.seed(33008)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- 1.5 * exp(0.6 * d$x) + 1.5 * (1 - exp(-0.6 * d$z))
  d$y <- rbinom(n, 1, plogis(lat - 2))
  d
}

test_that("the order of the terms does not change the fit", {
  d <- two_sat_data()
  f1 <- bf(y ~ b1 * (1 - exp(-b2 * z)) + b3 * (1 - exp(-b4 * x)),
           b1 + b2 + b3 + b4 ~ 1, nl = TRUE)
  f2 <- bf(y ~ b3 * (1 - exp(-b4 * x)) + b1 * (1 - exp(-b2 * z)),
           b3 + b4 + b1 + b2 ~ 1, nl = TRUE)
  a <- caught(frm(f1, family = bernoulli(), data = d))$value
  b <- caught(frm(f2, family = bernoulli(), data = d))$value
  la <- as.numeric(logLik(a))
  lb <- as.numeric(logLik(b))
  # a cap of 4 patterns by position read -163.5320 and -162.0191
  expect_equal(la, lb, tolerance = 1e-8)
})

test_that("a coefficient on a wide covariate starts at the body's scale", {
  # b2 near 1e-3 with xb in 0 to 3000 (the review's body bigx, seed 22);
  # starts of +-1 overflow exp(b2 * xb), and 0.79.1 stopped 19.0 short
  set.seed(22009)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- 0.5 * exp(6e-4 * d$xb)
  d$y <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf),
                     labels = 1:4))
  fm <- bf(y ~ b1 * exp(b2 * xb), b1 + b2 ~ 1, nl = TRUE)
  got <- caught(frm(fm, family = cumulative(), data = d))
  ref <- caught(frm(fm, family = cumulative(), data = d, start = list(
    beta = c("b1_(Intercept)" = -1, "b2_(Intercept)" = -1e-3))))$value
  ll <- as.numeric(logLik(got$value))
  ll_ref <- as.numeric(logLik(ref))
  # both stop at code 1 on this flat likelihood, a few 1e-3 apart
  expect_gte(ll, ll_ref - 1e-5 * abs(ll_ref))
  expect_true(has_class(got$conds, "frmtmb_nl_multistart"))
})

test_that("a prior's curvature enters a coefficient's unit", {
  # a * xs with xs about 1e-20 under normal(0, 1e-3): a unit read from the
  # likelihood alone stopped at code 1 below 0.79.1 (review, blocker 1)
  set.seed(6201)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$xs <- d$x * 1e-20
  body <- caught(frm(bf(y ~ a * xs, a ~ 1, nl = TRUE), data = d,
                     prior = set_prior("normal(0, 1e-3)", nlpar = "a")))
  lin <- frm(y ~ 0 + xs, data = d,
             prior = set_prior("normal(0, 1e-3)", class = "b"))
  expect_identical(body$value$opt$convergence, 0L)
  expect_false(any(vapply(body$conds, inherits, NA, "warning")))
  expect_equal(as.numeric(logLik(body$value)), as.numeric(logLik(lin)),
               tolerance = 1e-8)
})

test_that("the probe leaves a random-effect fit's objective where it was", {
  set.seed(7311)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), g = factor(sample(1:12, n, TRUE)))
  d$y <- (2 + rnorm(12, 0, 0.5)[d$g]) * exp(0.5 * d$x) + rnorm(n, 0, 0.3)
  fm <- bf(y ~ b1 * exp(b2 * x), b1 ~ 1 + (1 | g), b2 ~ 1, nl = TRUE)
  on <- caught(frm(fm, data = d))
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             caught(frm(fm, data = d)))
  expect_false(has_class(on$conds, "frmtmb_nl_multistart"))
  a <- on$value
  b <- off$value
  # the best point is what sdreport() and every later reader evaluate at;
  # without the restore it moved by up to 4.85e-7 here
  expect_identical(a$obj$env$last.par.best, b$obj$env$last.par.best)
  expect_identical(a$obj$env$last.par, b$obj$env$last.par)
  sa <- sdr_of(a)
  sb <- sdr_of(b)
  expect_identical(sa$par.fixed, sb$par.fixed)
  expect_identical(sa$cov.fixed, sb$cov.fixed)
})

test_that("a start the caller gives is not probed", {
  d <- sat_data(101)
  got <- caught(frm(sat_body, data = d, start = list(beta = c(
    "b1_(Intercept)" = 2, "b2_(Intercept)" = 0.5))))
  expect_false(has_class(got$conds, "frmtmb_nl_multistart"))
  expect_gt(got$value$estimates$beta[["b2_(Intercept)"]], 0)
})

test_that("the option turns the probe off, and the fit is then short", {
  d <- sat_data(101)
  on <- caught(frm(sat_body, data = d))$value
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             caught(frm(sat_body, data = d)))
  expect_false(has_class(off$conds, "frmtmb_nl_multistart"))
  ll_on <- as.numeric(logLik(on))
  # short by far more than any rounding: 71 log-likelihood units on 0.79.1
  expect_gt(ll_on - as.numeric(logLik(off$value)), 1e-3 * abs(ll_on))
})

test_that("a fit the probe cannot improve is the fit with the probe off", {
  set.seed(7302)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$y <- 2 * exp(0.5 * d$x) + rnorm(n, 0, 0.3)
  body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  # count the optimizations, so that the probe is seen to have run
  calls <- 0L
  orig <- optimize_obj
  local_mocked_bindings(optimize_obj = function(...) {
    calls <<- calls + 1L
    orig(...)
  })
  on <- caught(frm(body, data = d))
  n_on <- calls
  calls <- 0L
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             caught(frm(body, data = d)))
  expect_identical(n_on - calls, 4L)
  expect_false(has_class(on$conds, "frmtmb_nl_multistart"))
  expect_identical(on$value$opt$par, off$value$opt$par)
  expect_identical(as.numeric(logLik(on$value)),
                   as.numeric(logLik(off$value)))
  # the probes ran on the fit's own objective; the standard errors are
  # read from it afterwards and must be the ones at the returned optimum
  expect_equal(sqrt(diag(vcov(on$value))), sqrt(diag(vcov(off$value))),
               tolerance = 1e-8)
})

test_that("a body every coefficient enters as an added term is not probed", {
  set.seed(7303)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  body <- bf(y ~ a + b * x, a + b ~ 1, nl = TRUE)
  got <- caught(frm(body, data = d))
  expect_false(has_class(got$conds, "frmtmb_nl_multistart"))
  lin <- as.numeric(logLik(frm(y ~ x, data = d)))
  expect_equal(as.numeric(logLik(got$value)), lin, tolerance = 1e-8)
  frame <- got$value$frame
  lp <- frame[["linpreds"]][[which(vapply(frame[["linpreds"]], function(l) {
    !is.null(l[["nl_body"]])
  }, NA))]]
  bc <- nl_body_coefs(frame, lp)
  expect_true(all(nl_body_additive(frame, lp, got$value$estimates, bc$coefs)))
  # and a coefficient a product scales is not an added term
  body2 <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  f2 <- frm(body2, data = transform(d, y = exp(0.3 * x) + rnorm(n, 0, 0.2)))
  lp2 <- f2$frame[["linpreds"]][[which(vapply(f2$frame[["linpreds"]],
    function(l) !is.null(l[["nl_body"]]), NA))]]
  expect_false(any(nl_body_additive(f2$frame, lp2, f2$estimates,
                                    nl_body_coefs(f2$frame, lp2)$coefs)))
})

test_that("a coefficient on a covariate near 1e-20 is stepped in its units", {
  set.seed(6201)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$xt <- d$x * 1e-20
  lin <- frm(y ~ 0 + xt, data = d)
  body <- caught(frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d))
  ll_lin <- as.numeric(logLik(lin))
  # 0.79.1 stayed at the start, -299.33 against -272.50, at code 0
  expect_equal(as.numeric(logLik(body$value)), ll_lin, tolerance = 1e-8)
  expect_false(any(vapply(body$conds, inherits, NA, "warning")))
  expect_equal(body$value$estimates$beta[["a_(Intercept)"]],
               lin$estimates$beta[["xt"]], tolerance = 1e-6)
})

test_that("a body far from a large response reaches the linear maximum", {
  set.seed(9901)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
  lin <- as.numeric(logLik(frm(y ~ x, data = d)))
  body <- caught(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d))
  # 0.79.1 stopped at -3045.62 against -1658.76, at code 0
  expect_equal(as.numeric(logLik(body$value)), lin, tolerance = 1e-8)
  expect_false(any(vapply(body$conds, inherits, NA, "warning")))
})

test_that("a body of order-1 derivatives keeps its units at 1", {
  d <- sat_data(101)
  fit <- frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d)
  tpl <- fit$frame[["par_template"]]
  expect_null(nl_par_units(fit$frame, tpl, names(fit$obj$par)))
  d$xt <- d$x * 1e-20
  f2 <- frm(bf(y ~ a * xt, a ~ 1, nl = TRUE), data = d)
  u <- nl_par_units(f2$frame, f2$frame[["par_template"]], names(f2$obj$par))
  a_at <- which(names(f2$obj$par) == "beta")
  # one unit of `a` moves the body by about the response's own spread
  moved <- u[a_at] * sqrt(mean(d$xt^2)) / sd(d$y)
  expect_gt(moved, 0.5)
  expect_lt(moved, 2)
  expect_identical(log2(u[a_at]), round(log2(u[a_at])))
})

test_that("a prior whose curvature the unit does not read leaves no unit", {
  # gamma(2, 1000) on `a` in a * xs, xs about 1e-20: read as no prior, the
  # unit was 1.5e20 and the fit stopped at code 1, 992 below 0.79.1
  # (review, short check, R1)
  set.seed(6201)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$xs <- d$x * 1e-20
  got <- caught(frm(bf(y ~ a * xs, a ~ 1, nl = TRUE), data = d,
                    prior = set_prior("gamma(2, 1000)", nlpar = "a", lb = 0),
                    start = list(beta = c("a_(Intercept)" = 1))))
  fit <- got$value
  expect_identical(fit$opt$convergence, 0L)
  expect_false(any(vapply(got$conds, inherits, NA, "warning")))
  a_at <- which(names(fit$obj$par) == "beta")
  pu <- fit$par_units %||% rep(1, length(fit$obj$par))
  expect_identical(unname(pu[a_at]), 1)
})

test_that("a curve's coefficient is sized before the one that scales it", {
  # b1 * (1 - exp(-b2 * xs)) with xs about 1e-4 (review, short check, R2,
  # seed 1): b1 sized with b2 at 1 got b2's size, 32768, and 0.79.1 and
  # the first fix stayed 74 to 86 short with no warning
  set.seed(38)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$xs <- d$x * 1e-4
  d$y <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
  fm <- bf(y ~ b1 * (1 - exp(-b2 * xs)), b1 + b2 ~ 1, nl = TRUE)
  got <- caught(frm(fm, data = d))
  ref <- caught(frm(fm, data = d, start = list(beta = c(
    "b1_(Intercept)" = -1.5, "b2_(Intercept)" = -6000))))$value
  ll_ref <- as.numeric(logLik(ref))
  expect_gte(as.numeric(logLik(got$value)), ll_ref - 1e-8 * abs(ll_ref))
  expect_true(has_class(got$conds, "frmtmb_nl_multistart"))
})

test_that("a curve's coefficient on a covariate near 1e-6 fits as at 1", {
  # b2 near 1e6, stepped in absolute units, left this fit 45.9 short at
  # code 0 with no warning (review, short check 2, R3, poisson seed 3)
  set.seed(122)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), xp = round(runif(n, 0, 2), 3))
  d$y <- rpois(n, pmax(1.5 * exp(0.6 * d$x), 0.05))
  d$xs <- d$x * 1e-6
  at1 <- caught(frm(bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE),
                    family = poisson(), data = d))$value
  at6 <- caught(frm(bf(y ~ b1 * (1 - exp(-b2 * xs)), b1 + b2 ~ 1, nl = TRUE),
                    family = poisson(), data = d))
  l1 <- as.numeric(logLik(at1))
  expect_equal(as.numeric(logLik(at6$value)), l1, tolerance = 1e-6)
  expect_false(any(vapply(at6$conds, inherits, NA, "warning")))
  # one model with b2 rescaled, to a small part of b2's standard error
  # (under OpenBLAS the two estimates part by 1.7e-4 relative, at a
  # maximum the log-likelihood check above holds)
  fe <- fixef(at1)
  se_b2 <- fe[grep("^b2_", rownames(fe)), "Est.Error"]
  expect_lt(abs(at6$value$estimates$beta[["b2_(Intercept)"]] * 1e-6 -
                  at1$estimates$beta[["b2_(Intercept)"]]), 0.01 * se_b2)
})

test_that("a refit that did not settle is not returned silently", {
  d <- sat_data(101)
  # the refit's own record says it settled; here it is said not to, to
  # read what nl_fit_probed() then raises
  local_mocked_bindings(nl_sign_probe = function(run, frame, fit) {
    ll <- as.numeric(logLik(fit))
    list(fit = fit, held = list(), ll0 = ll - 1, ll1 = ll, n = 4L,
         settled = FALSE, ll_tape = ll + 1)
  })
  got <- caught(frm(sat_body, data = d))
  expect_true(has_class(got$conds, "frmtmb_nl_multistart_unsettled"))
  expect_true(has_class(got$conds, "frmtmb_nl_multistart"))
})

test_that("the sign probe runs without a basis whitened at the first stop", {
  # a Gompertz curve whose default fit stops at a lower maximum, where
  # nl_basis_descent() keeps a basis of the body's Jacobian; with that
  # basis the patterns' runs reached nothing (0.80.0 consolidation, 43.7
  # below the probe's point, review of lane nlstart's fresh row gomp 22)
  set.seed(22004)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  # the draws of the review's design, in its order
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  d$y <- 3 * exp(-2 * exp(-log(2) * d$xq)) + rnorm(n, 0, 0.5)
  fm <- bf(y ~ Asym * exp(-b2 * exp(-b3 * xq)), Asym ~ 1, b2 ~ 1, b3 ~ 1,
           nl = TRUE)
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             caught(frm(fm, data = d)))$value
  ob <- off$obj$env$.frmtmb_opt_basis
  # which point the first fit stops at depends on the optimizer's path;
  # the defect needs the body's basis on the tape
  skip_if_not(isTRUE(ob$nl), "the first fit kept no basis of the body")
  on <- caught(frm(fm, data = d))
  l_off <- as.numeric(logLik(off))
  l_on <- as.numeric(logLik(on$value))
  expect_gt(l_on - l_off, 1e-3 * abs(l_off))
  expect_true(has_class(on$conds, "frmtmb_nl_multistart"))
})

test_that("a probe winner on a ridge is not made silent by a second descent", {
  # binomial Michaelis-Menten (review of lane nlstart, fresh row mm 22):
  # the sign probe's winner runs to K -> Inf along the ridge to the linear
  # limit, 4.6 below the maximum, and its standard errors say so. A second
  # descent in the body's Jacobian basis from there slid further along
  # the ridge, the standard errors read finite and the fit was silent
  set.seed(22001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  # the draws of the review's design, in its order
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  d$nt <- 10L
  d$y <- rbinom(n, 10, plogis(3 * d$xp / (1 + d$xp) - 2))
  fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1, nl = TRUE)
  got <- caught(frm(fm, family = binomial(), data = d))
  at <- caught(frm(fm, family = binomial(), data = d, se = FALSE,
                   start = list(beta = c("Vm_(Intercept)" = -0.1,
                                         "K_(Intercept)" = -6))))
  l_def <- as.numeric(logLik(got$value))
  l_max <- as.numeric(logLik(at$value))
  warned <- any(vapply(got$conds, inherits, NA, "warning"))
  # the default fit reaches the maximum, or it says that it may not
  expect_true(l_max - l_def <= 1e-6 * abs(l_max) || warned)
})

test_that("the multistart message names the penalized value under a prior", {
  # merge step 3 of the 0.80.0 consolidation: what the fits maximized is
  # the penalized log-likelihood under a prior with a density, and the
  # message must say so and print it (release review of 0.80.0, item 11)
  fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  pr <- set_prior("normal(0, 5)", nlpar = "b1") +
    set_prior("normal(0, 5)", nlpar = "b2")
  set.seed(1)
  d <- data.frame(x = runif(200, 0, 4))
  d$y <- 3 * (1 - exp(-1.2 * d$x)) + rnorm(200, 0, 0.3)
  said <- function(got) {
    m <- Filter(function(cnd) {
      inherits(cnd, c("frmtmb_nl_multistart", "frmtmb_nl_restart"))
    }, got$conds)
    vapply(m, conditionMessage, "")
  }
  on <- caught(frm(fm, data = d, prior = pr))
  off <- caught(frm(fm, data = d))
  m_on <- said(on)
  m_off <- said(off)
  # whether a restart or the probe speaks depends on the optimizer's path
  skip_if(!length(m_on) || !length(m_off), "no restart or probe message")
  expect_true(all(grepl("penalized log-likelihood", m_on, fixed = TRUE)))
  expect_false(any(grepl("penalized", m_off, fixed = TRUE)))
  # the last value printed is the one the returned fit maximized
  nums <- regmatches(m_on, gregexpr("log-likelihood -?[0-9.]+", m_on))
  last <- as.numeric(sub(".*log-likelihood ", "", utils::tail(unlist(nums), 1)))
  lp <- as.numeric(logLik_penalized(on$value))
  expect_lt(abs(last / lp - 1), 1e-5)
})
