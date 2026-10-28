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
  }
  # each correlation has one slot, its class, as in brms, so no
  # internal route is named (class theta would refuse it)
  for (cl in c("cor", "cortime", "rescor")) {
    e <- tryCatch(set_prior("lkj(2)", class = cl, lb = 0),
                  error = conditionMessage)
    expect_match(e, "one prior slot", fixed = TRUE, info = cl)
    expect_no_match(e, "class = \"theta\"", fixed = TRUE, info = cl)
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
  # and names no route that is refused in turn: class theta on either
  # internal parameter is a second slot on the same coefficients
  e <- tryCatch(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = 0.1)), error = conditionMessage)
  expect_no_match(e, "class = \"theta\"", fixed = TRUE)
  expect_match(e, "order one only", fixed = TRUE)
  # the density still works at that order, which is what the Jacobian is
  # for
  expect_length(frmtmb:::resolve_prior_input(fit,
    set_prior("normal(0, 0.5)", class = "ar"))$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "theta", coef = "thetaac_2", lb = -1)),
    "internal parameter of ar", class = "frmtmb_error")
})

test_that("class theta refuses the covariance parameters a class reaches", {
  dd <- sim_ar()
  fm <- frm(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian(),
                 rescor = TRUE), data = dd)
  expect_error(frmtmb:::resolve_prior_input(fm,
    set_prior("", class = "theta", coef = "thetar_1", lb = -1)),
    "class = \"rescor\"", class = "frmtmb_error")
  fr <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(fr,
    set_prior("", class = "theta", coef = "theta_1", lb = -3)),
    "class = \"sd\", group = \"g\"", class = "frmtmb_error")
  # the class spelling bounds the standard deviation on its own scale
  re <- frmtmb:::resolve_prior_input(fr,
    set_prior("", class = "sd", lb = exp(-3)))
  expect_equal(unname(re$lower[["theta_1"]]), -3)

  # a name no covariance component carries lists the ones that exist
  f_ar <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(f_ar,
    set_prior("", class = "theta", coef = "nope_1", lb = 0)),
    "names no covariance parameter")
})

test_that("one residual-correlation parameter has one prior slot", {
  # through 0.75.0 default_prior() listed class ar and class theta coef
  # thetaac_1 for the same coefficient, and frm() applied whichever of
  # the two came last (dev/acprior-double.R)
  dd <- sim_ar(ng = 20)
  cases <- list(
    ar = list(y ~ x + ar(t, g), "ar"),
    arcov = list(y ~ x + ar(t, g, cov = TRUE), "ar"),
    ma = list(y ~ x + ma(t, g, cov = TRUE), "ma"),
    arma = list(y ~ x + arma(t, g, cov = TRUE), c("ar", "ma")),
    cosy = list(y ~ x + cosy(t, g), "cosy"),
    unstr = list(y ~ x + unstr(t, g), "cortime"))
  for (nm in names(cases)) {
    fo <- cases[[nm]][[1L]]
    cls <- cases[[nm]][[2L]]
    dp <- default_prior(fo, data = dd)
    expect_false(any(dp$class == "theta"), info = nm)
    expect_true(all(cls %in% dp$class), info = nm)
    fr <- frm(fo, data = dd, dry_run = "frame")
    pt <- fr[["par_template"]][["thetaac"]]
    for (k in seq_along(pt)) {
      e <- tryCatch(validate_prior(set_prior("normal(0, 1)",
                                             class = "theta",
                                             coef = paste0("thetaac_", k)),
                                   fo, data = dd),
                    error = identity)
      expect_s3_class(e, "frmtmb_error")
      # the refusal names the class that does reach the parameter
      own <- if (identical(nm, "arma")) cls[k] else cls
      expect_match(conditionMessage(e),
                   paste0("internal parameter of ", own, ","), fixed = TRUE,
                   info = paste(nm, k))
      expect_match(conditionMessage(e), paste0("class = \"", own, "\""),
                   fixed = TRUE, info = paste(nm, k))
    }
  }
  # frm() refuses it too, before any fitting
  expect_error(frm(y ~ x + cosy(t, g), data = dd,
                   prior = set_prior("", class = "theta",
                                     coef = "thetaac_1", lb = 0)),
               "internal parameter of cosy", class = "frmtmb_error")
  # the guard absent: the class itself takes the prior
  vp <- validate_prior(set_prior("normal(0.3, 0.2)", class = "ar"),
                       y ~ x + ar(t, g, cov = TRUE), data = dd)
  expect_true("normal(0.3, 0.2)" %in% vp$prior[vp$class == "ar"])
  # in a multivariate model the refusal names the response too
  e <- tryCatch(validate_prior(
    set_prior("normal(0, 1)", class = "theta", coef = "thetaac_2"),
    mvbf(bf(y ~ x + ar(t, g)), bf(y2 ~ x + ar(t, g))), data = dd),
    error = conditionMessage)
  expect_match(e, "class = \"ar\", resp = \"y2\"", fixed = TRUE)
})

test_that("the prior table bounds the classes as brms does", {
  # brms 2.23.0's default_prior() (dev/acprior-brms-defaults.R): lb -1
  # and ub 1 on a cov = TRUE ar and ma, none under cov = FALSE; cosy is
  # brms's [0, 1] on the sampling route only (frmtmb.sample's tests),
  # and frm()'s own window here
  dd <- sim_ar(ng = 20)
  bnd <- function(fo, cl) {
    dp <- default_prior(fo, data = dd)
    unlist(dp[dp$class == cl, c("lb", "ub")])
  }
  expect_equal(unname(bnd(y ~ x + ar(t, g, cov = TRUE), "ar")), c(-1, 1))
  expect_equal(unname(bnd(y ~ x + ma(t, g, cov = TRUE), "ma")), c(-1, 1))
  expect_true(all(is.na(bnd(y ~ x + ar(t, g), "ar"))))
  # an order-two coefficient is not on a box, and brms has no such term
  expect_true(all(is.na(bnd(y ~ x + ar(t, g, p = 2, cov = TRUE), "ar"))))
  # six time points: the window is (-1/5, 1)
  expect_equal(unname(bnd(y ~ x + cosy(t, g), "cosy")), c(-1 / 5, 1))
  # and the table passes back as prior =, its ends bounding nothing
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + cosy(t, g))) {
    f0 <- frm(fo, data = dd)
    f1 <- frm(fo, data = dd, prior = default_prior(fo, data = dd))
    expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  }
})

test_that("a bound at the end of a parameter's interval bounds nothing", {
  dd <- sim_ar()
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd)
  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = -1, ub = 1))
  # an end of the interval bounds nothing, so no bound is written
  expect_null(c(ri$lower, ri$upper))
  fc <- frm(bf(y ~ x + cosy(t, g)) + gaussian(), data = dd)
  rc <- frmtmb:::resolve_prior_input(fc,
    set_prior("", class = "cosy", lb = -1 / 5, ub = 1))
  expect_null(c(rc$lower, rc$upper))
  # the absent case: inside the interval it is a bound, and past its end
  # a refusal
  expect_true(is.finite(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = -0.5))$lower[[1L]]))
  expect_error(frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "ar", lb = -1.5)), "outside (-1, 1)",
    fixed = TRUE)
})

test_that("get_prior lists a row for every new class the model offers", {
  dd <- sim_ar()
  gp <- get_prior(bf(y ~ x + arma(t, g, cov = TRUE)) + gaussian(),
                  data = dd)
  expect_true(all(c("ar", "ma") %in% gp$class))
  # one row per parameter: the internal names of the ar and ma
  # coefficients are not a second spelling of them
  expect_false(any(gp$class == "theta"))

  gu <- get_prior(bf(y ~ x + unstr(t, g)) + gaussian(), data = dd)
  expect_true("cortime" %in% gu$class)
  expect_false(any(c("ar", "ma", "cosy") %in% gu$class))

  gm <- get_prior(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian(),
                       rescor = TRUE), data = dd)
  expect_true("rescor" %in% gm$class)
  # one row for the residual correlation, its class, as in brms
  expect_false(any(gm$class == "theta"))
})

test_that("prior_summary reports the new classes on a fit", {
  dd <- sim_ar()
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.5)", class = "ar"))
  txt <- paste(utils::capture.output(print(prior_summary(fit))),
               collapse = "\n")
  expect_match(txt, "ar ~ normal(0, 0.5)", fixed = TRUE)
})

test_that("a sigma prior under a cov = TRUE term is brms's, on sigma itself", {
  # brms's sigma under a cov = TRUE ar(), ma() or arma() is the
  # innovation sd (its cholesky_cor_*() are covariances), and so is
  # frmtmb's since sigma became the innovation sd. A typed
  # student_t(3, 0, 2.5) on sigma is brms's density at sigma = exp(eta)
  # with the log link's Jacobian, and no factor of the correlation
  # parameters enters it. frmtmb 0.76.0 held the marginal sd and read
  # the density at sigma * c with the Jacobian of c, so the penalty
  # moved with theta at a fixed eta.
  dd <- sim_ar(ng = 20)
  pr <- set_prior("student_t(3, 0, 2.5)", class = "sigma")
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                  y ~ x + arma(t, g, cov = TRUE),
                  y ~ x + ar(t, g, p = 2, cov = TRUE),
                  y ~ x + cosy(t, g), y ~ x + ar(t, g))) {
    fit <- frm(fo, data = dd)
    ac <- fit$frame$autocor[[1L]]
    ent <- frmtmb:::resolve_prior_input(fit, pr)$entries
    expect_length(ent, 1L)
    expect_null(ent[[1L]]$innov, label = deparse(fo))
    nlp <- frmtmb:::neg_log_prior_fn(ent)
    set.seed(5)
    v <- vapply(1:5, function(k) {
      p <- fit$frame$par_template
      p[["thetaac"]][] <- stats::rnorm(ac$npar)
      eta <- stats::rnorm(1, 0, 0.5)
      p[[ent[[1L]]$comp]][ent[[1L]]$idx] <- eta
      # brms's student_t_lpdf(sigma | 3, 0, 2.5), carried to eta
      ref <- stats::dt(exp(eta) / 2.5, 3, log = TRUE) - log(2.5) + eta
      c(nlp(p), ref)
    }, c(0, 0))
    expect_lt(diff(range(colSums(v))) / max(abs(v)), 1e-10,
              label = deparse(fo))
  }
})

test_that("a bound on sigma under a cov = TRUE term is brms's bound", {
  # brms takes lb/ub on its sigma, the innovation sd (stancode():
  # real<lower=0.2> sigma, dev/acprior-brms-sigmabound.R), which is the
  # parameter frmtmb holds, so the bound boxes it. frmtmb 0.76.0 refused
  # a lb above 0 or a ub there, because its sigma was the marginal sd.
  dd <- sim_ar(ng = 20)
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                  y ~ x + arma(t, g, cov = TRUE), y ~ x + cosy(t, g))) {
    fit <- frm(fo, data = dd)
    r <- frmtmb:::resolve_prior_input(fit,
      set_prior("normal(1, 1)", class = "sigma", lb = 0.2, ub = 3))
    expect_identical(unname(r$lower[["sigma_(Intercept)"]]), log(0.2),
                     label = deparse(fo))
    expect_identical(unname(r$upper[["sigma_(Intercept)"]]), log(3),
                     label = deparse(fo))
  }
  # a fit under the bound puts sigma, the innovation sd, inside it
  w <- "Standard errors are not available"
  fit <- allow_warnings(frm(y ~ x + ar(t, g, cov = TRUE), data = dd,
                            prior = set_prior("", class = "sigma", ub = 0.5)),
                        w, require = w)
  expect_lte(sigma(fit), 0.5 * (1 + 1e-8))
})
