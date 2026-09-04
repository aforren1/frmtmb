# The residual-correlation prior classes: brms's own class names for the
# surfaces frmtmb holds in `thetaac` and `thetar`.
#
# These parameters live on unconstrained internal scales, so a prior
# written about the NATURAL parameter (an AR coefficient, a cosy
# correlation) is carried inward with the log Jacobian of the map, the
# same change of variables class "sd" performs for a log standard
# deviation. The Jacobian is the part with no independent reference in
# the package, so it is checked against numeric differentiation first.

sim_ar <- function(seed = 3, ng = 60, k = 6, phi = 0.6) {
  set.seed(seed)
  n <- ng * k
  dd <- data.frame(x = stats::rnorm(n), t = rep(seq_len(k), ng),
                   g = factor(rep(seq_len(ng), each = k)))
  dd$y <- as.numeric(stats::arima.sim(list(ar = phi), n)) + 0.5 * dd$x
  dd$y2 <- stats::rnorm(n, 0.3 - 0.4 * dd$x, 1)
  dd
}

# the natural-scale estimate of a residual block, as summary() reports
# it: est_t is Fisher-z for anything living inside (-1, 1)
ac_natural <- function(fit) tanh(autocor_trans_rows(fit)$est_t)


test_that("the ARMA and cosy log Jacobians match numeric differentiation", {
  skip_if_not_installed("numDeriv")
  set.seed(9)
  # the map the density is carried through: unconstrained reals to the
  # AR (or MA) coefficients, via partial autocorrelations and Levinson
  f_lev <- function(t) autocor_levinson(autocor_pacf(t))
  for (p in 1:6) {
    for (rep in 1:3) {
      th <- stats::rnorm(p, 0, 1.2)
      num <- log(abs(det(numDeriv::jacobian(f_lev, th))))
      expect_equal(ac_trans_logjac(th, list(map = "levinson")), num,
                   tolerance = 1e-7, info = paste("order", p))
    }
  }
  for (d in c(3L, 5L)) {
    a <- 1 / (d - 1)
    th <- stats::rnorm(1)
    f_cosy <- function(t) -a + (1 + a) / (1 + exp(-t))
    num <- log(abs(numDeriv::grad(f_cosy, th)))
    expect_equal(ac_trans_logjac(th, list(map = "cosy", a = a)), num,
                 tolerance = 1e-7)
  }
})

test_that("class ar carries a density and a bound onto an ar() block", {
  dd <- sim_ar()
  form <- bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian()
  free <- frm(form, data = dd)
  phi0 <- ac_natural(free)[[1]]
  expect_lt(abs(phi0 - 0.6), 0.15)

  # the entry is one JOINT term over the block, because the map is not
  # elementwise
  e <- frmtmb:::resolve_prior_input(free,
    set_prior("normal(0, 0.5)", class = "ar"))$entries[[1L]]
  expect_identical(e$comp, "thetaac")
  expect_identical(e$dist$kind, "trans")

  # a bound bites exactly, on the natural scale the user wrote it on
  pinned <- suppressWarnings(
    frm(form, data = dd, prior = set_prior("", class = "ar", lb = 0.8)))
  expect_equal(ac_natural(pinned)[[1]], 0.8, tolerance = 1e-5)

  # and a tight density pulls the estimate toward its location
  shrunk <- frm(form, data = dd,
                prior = set_prior("normal(0, 0.05)", class = "ar"))
  expect_lt(ac_natural(shrunk)[[1]], phi0)
})

test_that("class cosy bounds map through the positive-definite window", {
  dd <- sim_ar()
  fit <- frm(bf(y ~ x + cosy(t, g)) + gaussian(), data = dd)
  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "cosy", lb = 0.1, ub = 0.5))
  expect_identical(names(ri$lower), "thetaac_1")
  # the map is monotone, so the internal bounds keep the order
  expect_lt(ri$lower[["thetaac_1"]], ri$upper[["thetaac_1"]])

  # a cosy correlation of d time points cannot go below -1/(d - 1), and
  # a bound outside that window is refused rather than clamped
  expect_error(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "cosy", lb = -0.9)), "outside the window")
})

test_that("class cortime and class rescor take lkj on their own matrix", {
  dd <- sim_ar()
  fu <- frm(bf(y ~ x + unstr(t, g)) + gaussian(), data = dd)
  eu <- frmtmb:::resolve_prior_input(fu,
    set_prior("lkj(2)", class = "cortime"))$entries[[1L]]
  expect_identical(eu$comp, "thetaac")
  expect_identical(eu$dist$kind, "lkj")
  # the whole time correlation, one entry per free correlation
  expect_length(eu$idx, 15L)

  fm <- frm(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian(),
                 rescor = TRUE), data = dd)
  em <- Filter(function(e) identical(e$comp, "thetar"),
               frmtmb:::resolve_prior_input(fm,
                 set_prior("lkj(4)", class = "rescor"))$entries)
  expect_length(em, 1L)
  expect_identical(em[[1L]]$dist$kind, "lkj")
  expect_identical(em[[1L]]$dist$d, 2L)
})

test_that("a class aimed at the wrong structure is refused by name", {
  dd <- sim_ar()
  f_ar <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(f_ar,
    set_prior("normal(0, 1)", class = "cosy")),
    "carries no \"cosy\" parameter")
  expect_error(frmtmb:::resolve_prior_input(f_ar,
    set_prior("normal(0, 1)", class = "ma")), "carries no \"ma\"")
  # a model with no residual structure at all says so differently
  plain <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(plain,
    set_prior("normal(0, 1)", class = "ar")),
    "no residual autocorrelation term")
  expect_error(frmtmb:::resolve_prior_input(plain,
    set_prior("lkj(2)", class = "rescor")), "This model has none")
})

test_that("the matrix-valued classes refuse lb/ub, naming the alternative", {
  for (cl in c("cor", "cortime", "rescor")) {
    expect_error(set_prior("lkj(2)", class = cl, lb = 0),
                 "takes no lb/ub", info = cl)
    expect_error(set_prior("lkj(2)", class = cl, lb = 0),
                 "thetaac_1", info = cl)
  }
  # and each of them still demands the density that fits it
  expect_error(set_prior("normal(0, 1)", class = "rescor"),
               "takes an lkj\\(\\) prior")
  expect_error(set_prior("lkj(2)", class = "ar"),
               "belongs to class")
})

test_that("a higher-order ar coefficient refuses a bound, with the reason", {
  dd <- sim_ar()
  fit <- frm(bf(y ~ x + ar(t, g, p = 2, cov = TRUE)) + gaussian(),
             data = dd)
  # ar[1] is a function of both internal parameters, so no box in
  # internal space is the box the user asked for
  expect_error(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = 0.1)), "takes no lb/ub at order 2")
  expect_error(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = 0.1)), "class = \"theta\"")
  # the density still works at that order, which is what the Jacobian is
  # for, and the internal escape hatch still bounds one parameter
  expect_length(frmtmb:::resolve_prior_input(fit,
    set_prior("normal(0, 0.5)", class = "ar"))$entries, 1L)
  expect_identical(names(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "theta", coef = "thetaac_2", lb = -1))$lower),
    "thetaac_2")
})

test_that("class theta reaches every covariance component by name", {
  dd <- sim_ar()
  f_ar <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd)
  ri <- frmtmb:::resolve_prior_input(f_ar,
    set_prior("", class = "theta", coef = "thetaac_1", lb = -2, ub = 2))
  expect_identical(names(ri$lower), "thetaac_1")
  expect_identical(unname(ri$lower), -2)

  fm <- frm(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian(),
                 rescor = TRUE), data = dd)
  rr <- frmtmb:::resolve_prior_input(fm,
    set_prior("", class = "theta", coef = "thetar_1", lb = -1))
  expect_identical(names(rr$lower), "thetar_1")

  # a name no covariance component carries lists the ones that exist
  expect_error(frmtmb:::resolve_prior_input(f_ar,
    set_prior("", class = "theta", coef = "nope_1", lb = 0)),
    "names no covariance parameter")
})

test_that("get_prior lists a row for every new class the model offers", {
  dd <- sim_ar()
  gp <- get_prior(bf(y ~ x + arma(t, g, cov = TRUE)) + gaussian(),
                  data = dd)
  expect_true(all(c("ar", "ma") %in% gp$class))
  # per-element internal rows, so every bound has a spelling to copy
  expect_true(all(c("thetaac_1", "thetaac_2") %in%
                    gp$coef[gp$class == "theta"]))
  # a model with no random effects offers no bare theta row: `$theta`
  # partial-matches `thetaac`, which once listed one
  expect_false(any(gp$class == "theta" & !nzchar(gp$coef)))

  gu <- get_prior(bf(y ~ x + unstr(t, g)) + gaussian(), data = dd)
  expect_true("cortime" %in% gu$class)
  expect_false(any(c("ar", "ma", "cosy") %in% gu$class))

  gm <- get_prior(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian(),
                       rescor = TRUE), data = dd)
  expect_true("rescor" %in% gm$class)
  expect_true("thetar_1" %in% gm$coef[gm$class == "theta"])
})

test_that("prior_summary reports the new classes on a fit", {
  dd <- sim_ar()
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.5)", class = "ar"))
  txt <- paste(utils::capture.output(print(prior_summary(fit))),
               collapse = "\n")
  expect_match(txt, "class=ar")
  expect_match(txt, "normal(0, 0.5)", fixed = TRUE)
})
