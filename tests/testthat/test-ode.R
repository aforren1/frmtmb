# frm_ode(): the per-group ODE solve inside a nonlinear predictor.
#
# The reference numbers come from dev/ode/probeA2-frmtmb-nl.R, which
# fits the same 12-subject, 8-timepoint population PK data with a
# hand-written per-subject solve loop, and which was itself cross-
# checked against a hand-rolled RTMB MakeADFun and against nlmixr2
# FOCEi (dev/ode-feasibility.md, section 2).
#
# The order of the tests in this file matters. A deliberately failing
# ODE solve can poison later MakeADFun objects in the same session
# (dev/ode-feasibility.md, section 7), so every test that fails a solve
# on purpose lives at the end.

# One-compartment first-order oral absorption, single bolus into the
# depot:  dA/dt = -ka A,  dC/dt = ka A / V - ke C.
pk_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

pk_analytic <- function(t, ka, ke, V, D) {
  D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
}

# The probe's simulator, verbatim in effect: seed 2026, 12 subjects,
# 8 timepoints, sd(lka) = 0.30, sd(lke) = 0.25, sigma = 0.30.
sim_pk <- function(n_id = 12, seed = 2026) {
  set.seed(seed)
  tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  id <- factor(rep(seq_len(n_id), each = length(tt)))
  time <- rep(tt, n_id)
  b_ka <- stats::rnorm(n_id, 0, 0.30)
  b_ke <- stats::rnorm(n_id, 0, 0.25)
  ka <- exp(0 + b_ka)[as.integer(id)]
  ke <- exp(log(0.2) + b_ke)[as.integer(id)]
  mu <- pk_analytic(time, ka, ke, 10, 100)
  data.frame(id = id, time = time, dose = 100,
             conc = mu + stats::rnorm(length(mu), 0, 0.30))
}

pk_form <- bf(
  conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                 parms = list(exp(lka), exp(lke), exp(lV)), group = id,
                 states = c("depot", "central"), output = "central"),
  lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE
)

# --- checks that need no solver ------------------------------------

test_that("a function argument of an nl body is not looked for in data", {
  # frm_ode(pk_dyn, ...) puts `pk_dyn` in all.vars() of the body; asking
  # model.frame() for a column of that name is the wrong question
  d <- sim_pk(3)
  fr <- frm(pk_form + gaussian(), data = d, dry_run = "frame",
            start = list(beta = c(0, log(0.25), log(8))))
  nl <- fr$linpreds[["conc.mu"]]
  expect_false("pk_dyn" %in% names(nl$data_list))
  expect_setequal(names(nl$data_list), c("dose", "time", "id"))
})

test_that("a misspelled column that matches a function is named", {
  # the cost of resolving body names outside the data: a typo whose
  # name happens to be a base function is no longer caught by
  # model.frame(), so the failure has to be re-raised with the suspects
  d <- sim_pk(3)
  names(d)[names(d) == "time"] <- "t"     # the column the body wants
  form <- bf(conc ~ b0 * t, b0 ~ 1, nl = TRUE)   # ... but data has `t`
  expect_s3_class(frm(form + gaussian(), data = d, dry_run = "frame"),
                  "frmtmb_frame")

  names(d)[names(d) == "t"] <- "tim"      # now `t` is base::t, a typo
  expect_error(frm(form + gaussian(), data = d, start = list(beta = 1)),
               "resolved outside the data")
  expect_error(frm(form + gaussian(), data = d, start = list(beta = 1)),
               "\\bt\\b")
})

test_that("a column wins over a same-named function", {
  d <- sim_pk(3)
  d$c <- d$dose            # shadows base::c
  form <- bf(conc ~ frm_ode(pk_dyn, init = list(c, 0), times = time,
                            parms = list(exp(lka), exp(lke), exp(lV)),
                            group = id, output = 2L),
             lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE)
  fr <- frm(form + gaussian(), data = d, dry_run = "frame",
            start = list(beta = c(0, log(0.25), log(8))))
  expect_true("c" %in% names(fr$linpreds[["conc.mu"]]$data_list))
})

test_that("a dynamics parameter that varies within a group is refused", {
  # the trap of dev/ode-feasibility.md section 5: without this the fit
  # leaves the coefficient at its start value with an indefinite Hessian
  d <- sim_pk(4)
  d$phase <- factor(ifelse(d$time <= 2, "early", "late"))
  form <- bf(
    conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                   parms = list(exp(lka), exp(lke), exp(lV)), group = id,
                   output = 2L),
    lka ~ 1 + (1 | id), lke ~ 1 + phase + (1 | id), lV ~ 1, nl = TRUE
  )
  expect_error(
    frm(form + gaussian(), data = d, dry_run = "frame",
        start = list(beta = c(0, log(0.25), 0, log(8)))),
    "not constant within 'id'"
  )
})

test_that("a covariate constant within group is allowed", {
  d <- sim_pk(4)
  set.seed(1)
  d$wt <- rep(stats::rnorm(4), each = 8)
  form <- bf(
    conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                   parms = list(exp(lka), exp(lke), exp(lV)), group = id,
                   output = 2L),
    lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1 + wt, nl = TRUE
  )
  expect_s3_class(
    frm(form + gaussian(), data = d, dry_run = "frame",
        start = list(beta = c(0, log(0.25), log(8), 0))),
    "frmtmb_frame"
  )
})

test_that("a missing backend names the r-universe repository", {
  # this is the CRAN case: RTMBode is not on CRAN and the check farm
  # will not install it, so the error has to say where it comes from
  local_mocked_bindings(ode_has_pkg = function(pkg) pkg != "RTMBode")
  expect_error(frm_ode(function(t, y, p) list(-y), init = list(1),
                       times = 1, parms = list(1)),
               "kaskr\\.r-universe\\.dev")
  local_mocked_bindings(ode_has_pkg = function(pkg) FALSE)
  expect_error(frm_ode(function(t, y, p) list(-y), init = list(1),
                       times = 1, parms = list(1)),
               "RTMBode and deSolve packages")
})

# --- checks that need RTMBode ---------------------------------------

test_that("the solver reproduces the closed-form solution", {
  skip_if_not_installed("RTMBode")
  d <- sim_pk()
  n <- nrow(d)
  got <- frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
                 parms = list(rep(1, n), rep(0.2, n), rep(10, n)),
                 group = d$id, states = c("depot", "central"),
                 output = "central")
  expect_equal(got, pk_analytic(d$time, 1, 0.2, 10, 100),
               tolerance = 1e-6)
})

test_that("the default returns every state as a matrix", {
  skip_if_not_installed("RTMBode")
  d <- sim_pk(3)
  m <- frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
               parms = list(1, 0.2, 10), group = d$id)
  expect_identical(dim(m), c(nrow(d), 2L))
  expect_equal(m[, 1], 100 * exp(-d$time), tolerance = 1e-6)
  expect_equal(m[, 2], pk_analytic(d$time, 1, 0.2, 10, 100),
               tolerance = 1e-6)
})

test_that("rows are scattered back in input order", {
  skip_if_not_installed("RTMBode")
  d <- sim_pk(4)
  set.seed(3)
  sh <- d[sample(nrow(d)), ]
  got <- frm_ode(pk_dyn, init = list(sh$dose, 0), times = sh$time,
                 parms = list(1, 0.2, 10), group = sh$id, output = 2L)
  expect_equal(got, pk_analytic(sh$time, 1, 0.2, 10, 100),
               tolerance = 1e-6)
})

test_that("duplicate times, an observation at t0, and ragged groups work", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(
    id = factor(c(1, 1, 1, 1, 2, 2, 2)),
    time = c(0, 1, 1, 4, 0.5, 2, 8),
    dose = 100
  )
  got <- frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id, output = 2L)
  expect_equal(got[1], 0)                       # C(0) = 0
  expect_equal(got[2], got[3])                  # the duplicate time
  expect_equal(got, pk_analytic(d$time, 1, 0.2, 10, 100),
               tolerance = 1e-6)
})

test_that("t0 shifts the origin", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(1), time = c(1.25, 3), dose = 100)
  got <- frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id, t0 = 1,
                 output = 2L)
  expect_equal(got, pk_analytic(d$time - 1, 1, 0.2, 10, 100),
               tolerance = 1e-6)
  expect_error(
    frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id, t0 = 2),
    "before t0"
  )
})

test_that("group = NULL solves the whole data as one system", {
  skip_if_not_installed("RTMBode")
  tt <- c(0.25, 1, 4, 12)
  got <- frm_ode(pk_dyn, init = list(100, 0), times = tt,
                 parms = list(1, 0.2, 10), output = 2L)
  expect_equal(got, pk_analytic(tt, 1, 0.2, 10, 100), tolerance = 1e-6)
})

test_that("the arguments are validated", {
  skip_if_not_installed("RTMBode")
  d <- sim_pk(2)
  expect_error(
    frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id, method = "rk4"),
    "adaptive integrator"
  )
  expect_error(
    frm_ode(pk_dyn, init = c(100, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id),
    "pass a list"
  )
  expect_error(
    frm_ode(pk_dyn, init = list(d$time, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id),
    "not constant within group"
  )
  expect_error(
    frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id,
            states = c("a", "b", "c")),
    "`states` names 3 states"
  )
  expect_error(
    frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id,
            states = c("a", "b"), output = "central"),
    "not in `states`"
  )
  expect_error(
    frm_ode(pk_dyn, init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id, output = 5L),
    "index states 1 to 2"
  )
  expect_error(
    frm_ode("not a function", init = list(d$dose, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id),
    "must be a function"
  )
})

test_that("a large single system warns about the Laplace ceiling", {
  skip_if_not_installed("RTMBode")
  decay <- function(t, y, p) list(-p[1] * y)
  expect_warning(
    frm_ode(decay, init = as.list(rep(1, 8)), times = c(0.5, 1),
            parms = list(0.3), output = 1L),
    "Laplace"
  )
})

test_that("the population PK fit matches the probe reference", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  d <- sim_pk()
  fit <- frm(pk_form + gaussian(), data = d,
             start = list(beta = c(0, log(0.25), log(8))), se = TRUE)

  # dev/ode/probeA2-frmtmb-nl.R, cross-checked against a hand-rolled
  # RTMB MakeADFun (same objective to 6 digits) and nlmixr2 FOCEi
  expect_equal(as.numeric(logLik(fit)), -60.462931, tolerance = 1e-6)
  expect_equal(unname(fixef(fit)$lka), -0.2292, tolerance = 1e-3)
  expect_equal(unname(fixef(fit)$lke), -1.5932, tolerance = 1e-3)
  expect_equal(unname(fixef(fit)$lV), 2.2704, tolerance = 1e-3)
  expect_equal(sigma(fit), 0.3040, tolerance = 1e-3)
  expect_equal(AIC(fit), 132.926, tolerance = 1e-5)

  sds <- vapply(VarCorr(fit), function(v) sqrt(v[1L, 1L]), 0)
  expect_equal(unname(sds), c(0.2533, 0.2624), tolerance = 1e-3)

  # standard errors exist and are finite: the Hessian is positive
  # definite, which the within-group-covariate trap destroys
  expect_true(all(is.finite(summary(fit)$coefficients$lka[, "Std. Error"])))

  # the random effects recover the simulated subject deviations
  set.seed(2026)
  b_ka <- stats::rnorm(12, 0, 0.30)
  expect_gt(stats::cor(as.numeric(ranef(fit)[[1L]]), b_ka), 0.95)
})

test_that("predict() re-solves on newdata", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  d <- sim_pk(6)
  fit <- frm(pk_form + gaussian(), data = d,
             start = list(beta = c(0, log(0.25), log(8))))
  p_in <- predict(fit)
  keep <- d$id %in% c("1", "3")
  expect_equal(predict(fit, newdata = d[keep, ]), unname(p_in[keep]),
               tolerance = 1e-6, ignore_attr = TRUE)
  # a denser grid the fit never saw
  nd <- expand.grid(time = seq(0.5, 12, by = 0.5),
                    id = factor("1", levels = levels(d$id)))
  nd$dose <- 100
  expect_length(predict(fit, newdata = nd), nrow(nd))
})

test_that("an integrator choice does not move the likelihood", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  d <- sim_pk(6)
  ll <- vapply(c("lsoda", "adams", "ode45"), function(m) {
    # the method has to be a literal in the body: every symbol of a
    # nonlinear body that is not a nonlinear parameter is looked for as
    # a data column
    form <- eval(bquote(bf(
      conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                     parms = list(exp(lka), exp(lke), exp(lV)),
                     group = id, output = 2L, method = .(m)),
      lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE
    )))
    # convergence warnings differ across BLAS builds (mac reports false
    # convergence at gradient 2e-3); the logLik agreement below is the
    # assertion that matters
    as.numeric(logLik(suppressWarnings(
      frm(form + gaussian(), data = d,
          start = list(beta = c(0, log(0.25), log(8))))
    )))
  }, 0)
  expect_equal(unname(diff(range(ll))), 0, tolerance = 1e-5)
})

# --- one solve, a different state per row ---------------------------

# Lotka-Volterra: two species observed on the same response column, one
# row per species per time. dev/ode-feasibility.md section 10.
lv_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(p[1] * y[1] - p[2] * y[1] * y[2],
         p[3] * y[1] * y[2] - p[4] * y[2]))
}

test_that("a per-row output is the two-call spelling, exactly", {
  skip_if_not_installed("RTMBode")
  tt <- rep(c(1, 3, 5, 8, 12), each = 2)
  cmt <- rep(c(1L, 2L), 5)
  pv <- list(0.55, 0.028, 0.026, 0.84)
  one <- frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
                 states = c("hare", "lynx"), output = cmt,
                 atol = 1e-10, rtol = 1e-10)
  a <- frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
               states = c("hare", "lynx"), output = "hare",
               atol = 1e-10, rtol = 1e-10)
  b <- frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
               states = c("hare", "lynx"), output = "lynx",
               atol = 1e-10, rtol = 1e-10)
  two <- ifelse(cmt == 1L, a, b)
  expect_identical(one, two)

  # by name as well as by position
  by_name <- frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
                     states = c("hare", "lynx"),
                     output = c("hare", "lynx")[cmt],
                     atol = 1e-10, rtol = 1e-10)
  expect_identical(by_name, one)

  # a selection of states is still a selection of states: `output` is
  # read per row only when there are more rows than states
  m <- frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
               output = c(1L, 2L), atol = 1e-10, rtol = 1e-10)
  expect_identical(dim(m), c(length(tt), 2L))

  expect_error(
    frm_ode(lv_dyn, init = list(30, 4), times = tt, parms = pv,
            output = rep(3L, length(tt))),
    "index states 1 to 2"
  )
})

test_that("a per-row output shares one solve across groups", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(rep(c("a", "b"), each = 6)),
                  time = rep(rep(c(1, 4, 9), each = 2), 2),
                  cmt = rep(c(1L, 2L), 6))
  d$k <- rep(c(0.3, 0.6), each = 6)
  one <- frm_ode(pk_dyn, init = list(100, 0), times = d$time,
                 parms = list(1, d$k, 10), group = d$id, output = d$cmt,
                 atol = 1e-10, rtol = 1e-10)
  m <- frm_ode(pk_dyn, init = list(100, 0), times = d$time,
               parms = list(1, d$k, 10), group = d$id,
               atol = 1e-10, rtol = 1e-10)
  expect_identical(one, ifelse(d$cmt == 1L, m[, 1], m[, 2]))
})

test_that("a parent-and-metabolite fit is the same one call or two", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # depot -> parent -> metabolite, both measured on one assay column
  pm_dyn <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[1] * y[1],
           p[1] * y[1] - (p[2] + p[3]) * y[2],
           p[2] * y[2] - p[4] * y[3]))
  }
  set.seed(9)
  tt <- c(0.5, 1, 2, 4, 8, 12)
  n_id <- 5
  base <- expand.grid(time = tt, cmt = c(2L, 3L))
  d <- do.call(rbind, lapply(seq_len(n_id), function(i)
    data.frame(id = factor(i, levels = seq_len(n_id)), base)))
  truth <- frm_ode(pm_dyn, init = list(100, 0, 0), times = d$time,
                   parms = list(1.1, 0.35, 0.2, 0.25), group = d$id,
                   output = d$cmt, atol = 1e-10, rtol = 1e-10)
  d$y <- truth + stats::rnorm(nrow(d), 0, 1.5)
  d$is_p <- as.numeric(d$cmt == 2L)
  st <- list(beta = c(0, log(0.3), log(0.2), log(0.3)))

  one <- frm(bf(y ~ frm_ode(pm_dyn, init = list(100, 0, 0), times = time,
                            parms = list(exp(lka), exp(lkm), exp(lke),
                                         exp(lkme)),
                            group = id, output = cmt),
                lka ~ 1, lkm ~ 1, lke ~ 1, lkme ~ 1, nl = TRUE) +
               gaussian(), data = d, start = st)
  two <- frm(bf(y ~ is_p * frm_ode(pm_dyn, init = list(100, 0, 0),
                                   times = time,
                                   parms = list(exp(lka), exp(lkm),
                                                exp(lke), exp(lkme)),
                                   group = id, output = 2L) +
                  (1 - is_p) * frm_ode(pm_dyn, init = list(100, 0, 0),
                                       times = time,
                                       parms = list(exp(lka), exp(lkm),
                                                    exp(lke), exp(lkme)),
                                       group = id, output = 3L),
                lka ~ 1, lkm ~ 1, lke ~ 1, lkme ~ 1, nl = TRUE) +
               gaussian(), data = d, start = st)
  expect_equal(as.numeric(logLik(one)), as.numeric(logLik(two)),
               tolerance = 1e-8)
  expect_equal(unlist(fixef(one)), unlist(fixef(two)), tolerance = 1e-6)
  # and it recovers what generated the data
  fx <- unlist(fixef(one))
  expect_equal(unname(exp(fx[["lka.(Intercept)"]])), 1.1, tolerance = 0.3)
  expect_equal(unname(exp(fx[["lkme.(Intercept)"]])), 0.25,
               tolerance = 0.3)
})

# --- a combined proportional-plus-additive error model ---------------

test_that("nlf(sigma ~ ) writes the combined error model", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # sd = sqrt(add^2 + (prop * mu)^2), which on the log link sigma is
  # reported on is 0.5 * log(exp(2 la) + exp(2 lp) mu^2)
  set.seed(13)
  n_id <- 10
  tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  d <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                  time = rep(tt, n_id), dose = 100)
  ka <- exp(stats::rnorm(n_id, 0, 0.2))[as.integer(d$id)]
  ke <- exp(stats::rnorm(n_id, log(0.2), 0.2))[as.integer(d$id)]
  mu <- pk_analytic(d$time, ka, ke, 10, 100)
  add <- 0.15
  prop <- 0.10
  d$conc <- mu + stats::rnorm(nrow(d), 0,
                              sqrt(add^2 + (prop * mu)^2))

  fit <- frm(
    bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                      parms = list(exp(lka), exp(lke), exp(lV)),
                      group = id, output = 2L),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      nlf(sigma ~ 0.5 * log(exp(2 * ladd) + exp(2 * lprop) * mu^2)) +
      lf(ladd ~ 1, lprop ~ 1) + gaussian(),
    data = d, start = list(beta = c(0, log(0.25), log(8)),
                           betad = c(log(0.2), log(0.2))))

  fx <- unlist(fixef(fit))
  expect_equal(unname(exp(fx[["ladd.(Intercept)"]])), add, tolerance = 0.8)
  expect_equal(unname(exp(fx[["lprop.(Intercept)"]])), prop,
               tolerance = 0.8)
  # the residual sd rises with the fitted value, which a single sigma
  # cannot do
  s <- predict(fit, dpar = "sigma")
  p <- predict(fit)
  expect_gt(stats::cor(as.numeric(s), as.numeric(p)), 0.9)

  flat <- frm(
    bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                      parms = list(exp(lka), exp(lke), exp(lV)),
                      group = id, output = 2L),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = d, start = list(beta = c(0, log(0.25), log(8))))
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(flat)))
})

# --- deliberately failing solves, last in the file ------------------

test_that("a failed solve becomes a finite penalty or an error", {
  skip_if_not_installed("RTMBode")
  # dy/dt = p y^2 reaches infinity at t = 1/(p y0). Group 1 is asked for
  # times well inside that; group 2 is asked for times past it, which is
  # what the optimizer does to a real model when it probes an extreme
  # rate constant.
  blowup <- function(t, y, p) list(p[1] * y * y)
  d <- data.frame(id = factor(c(1, 1, 2, 2)), time = c(1, 2, 1, 2))
  parms <- list(c(0.01, 0.01, 20, 20))

  got <- suppressWarnings(
    frm_ode(blowup, init = list(1), times = d$time, parms = parms,
            group = d$id, output = 1L, method = "ode45"))
  expect_equal(got[1:2], 1 / (1 - 0.01 * c(1, 2)), tolerance = 1e-5)
  expect_equal(got[3:4], c(1e6, 1e6))
  expect_true(all(is.finite(got)))

  got2 <- suppressWarnings(
    frm_ode(blowup, init = list(1), times = d$time, parms = parms,
            group = d$id, output = 1L, method = "ode45", penalty = 42))
  expect_equal(got2[3:4], c(42, 42))

  expect_error(
    suppressWarnings(
      frm_ode(blowup, init = list(1), times = d$time, parms = parms,
              group = d$id, output = 1L, method = "ode45", on_error = "error")),
    "failed to solve group '2'"
  )
})

test_that("a penalty is never written silently", {
  skip_if_not_installed("RTMBode")
  blowup <- function(t, y, p) list(p[1] * y * y)
  d <- data.frame(id = factor(c(1, 1, 2, 2)), time = c(1, 2, 1, 2))
  parms <- list(c(0.01, 0.01, 20, 20))

  expect_warning(
    got <- frm_ode(blowup, init = list(1), times = d$time, parms = parms,
                   group = d$id, output = 1L, method = "ode45"),
    "the solve failed for 1 of 2 groups \\(2\\)"
  )
  expect_equal(got[3:4], c(1e6, 1e6))

  # and it can be read back after the warning was suppressed
  log <- frm_ode_failures()
  expect_identical(log$groups, "2")
  expect_identical(log$n_groups, 2L)
  expect_identical(log$penalty, 1e6)

  # a clean call clears the record
  frm_ode(blowup, init = list(1), times = d$time[1:2],
          parms = list(0.01), group = d$id[1:2], output = 1L)
  expect_null(frm_ode_failures())
})

test_that("a numeric-only body is loud when a solve fails at tape time", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # init and parms hold no estimated parameter, so the body is evaluated
  # numerically as the tape is built. The penalty would otherwise be
  # baked into the objective with nothing said.
  blowup <- function(t, y, p) list(p[1] * y * y)
  set.seed(4)
  d <- data.frame(id = factor(rep(c(1, 2, 3), each = 3)),
                  time = rep(c(1, 2, 3), 3),
                  p = rep(c(0.01, 0.02, 20), each = 3))
  d$conc <- 1 / (1 - d$p * d$time) + stats::rnorm(9, 0, 0.01)

  form <- bf(conc ~ b0 * frm_ode(blowup, init = list(1), times = time,
                                 parms = list(p), group = id, output = 1L),
             b0 ~ 1, nl = TRUE)
  # capture_warnings: deSolve's DLSODA give-up warnings propagate by
  # design and would otherwise land in the test record
  ws <- capture_warnings(frm(form + gaussian(), data = d,
                             start = list(beta = 1)))
  expect_true(any(grepl("the solve failed for 1 of 3 groups \\(3\\)", ws)))
  expect_identical(frm_ode_failures()$groups, "3")
})

test_that("a penalty reached through predict() warns and names the group", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  blowup <- function(t, y, p) list(p[1] * y * y)
  set.seed(5)
  d <- data.frame(id = factor(rep(c(1, 2, 3), each = 3)),
                  time = rep(c(1, 2, 3), 3))
  d$conc <- 1 / (1 - 0.05 * d$time) + stats::rnorm(9, 0, 0.01)

  form <- bf(conc ~ frm_ode(blowup, init = list(1), times = time,
                            parms = list(exp(lp)), group = id,
                            output = 1L),
             lp ~ 1, nl = TRUE)
  fit <- frm(form + gaussian(), data = d, start = list(beta = log(0.05)))

  # a time past the singularity at t = 1 / p, which the fit never saw
  nd <- data.frame(id = factor(c("1", "2"), levels = levels(d$id)),
                   time = c(2, 500))
  # capture_warnings absorbs the propagated DLSODA give-up warnings too
  ws <- capture_warnings(p <- predict(fit, newdata = nd))
  expect_true(any(grepl("the solve failed for 1 of 2 groups \\(2\\)", ws)))
  expect_equal(unname(p[2]), 1e6)
  expect_identical(frm_ode_failures()$groups, "2")
})

test_that("dynamics need no ADoverload boilerplate of their own", {
  skip_if_not_installed("RTMBode")
  # base c() and dy[i] <- in a user function, no ADoverload lines: the
  # wrap in frm_ode() puts the overloads in scope, and a pre-compiled
  # function is the worst case for it
  set.seed(41)
  d1 <- data.frame(t = rep(seq(0.5, 6, by = 0.5), 4),
                   id = factor(rep(1:4, each = 12)))
  d1$y <- 10 * exp(-0.3 * d1$t) * exp(stats::rnorm(nrow(d1), 0, 0.05))
  bare <- compiler::cmpfun(function(t, y, p) {
    dy <- numeric(1)
    dy[1] <- -p[1] * y[1]
    list(c(dy[1]))
  })
  boiler <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    dy <- numeric(1)
    dy[1] <- -p[1] * y[1]
    list(c(dy[1]))
  }
  mk <- function(dyn) frm(
    bf(y ~ log(frm_ode(dyn, init = list(exp(lv0)), times = t,
                       parms = list(exp(lk)), group = id)),
       lk ~ 1, lv0 ~ 1, nl = TRUE),
    family = lognormal(), data = d1,
    start = list(beta = log(c(0.3, 10))))
  f1 <- mk(bare)
  f2 <- mk(boiler)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
  expect_equal(unname(unlist(fixef(f1))), unname(unlist(fixef(f2))),
               tolerance = 1e-8)
})
