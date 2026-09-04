# frm_ode(tv = ): piecewise-constant time-varying inputs.
#
# This is the one place where a dynamics input MAY vary inside a solve
# group. It rides the segmented-solve machinery that dosing events
# already use: the solve is split at each within-group change point, the
# breakpoints are unioned with the event times, and each segment's
# dynamics see that segment's value as an extra parameter. Parameters
# are tape inputs, so an estimated time-varying value (a covariate-
# dependent clearance that changes at a known time) is differentiated
# exactly, while the CHANGE POINTS stay data - they decide where the
# solve is split, which is settled before the tape is built.
#
# Semantics: last observation carried forward, which is what rxode2 does
# with covariates. A row's value is in force from that row's time until
# the next change. The state itself is continuous across a change point
# (a covariate moves the derivative, not the state), so an observation
# exactly at a change point reads the value the PRE-change dynamics
# produced. That is pinned below.
#
# As in the other ODE files, anything that fails a solve on purpose
# lives at the end and captures its warnings.

tv_decay <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}

# one-compartment oral absorption whose elimination rate is the
# time-varying input: parms = (ka, V), tv = (ke)
tv_pk <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] / p[2] - p[3] * y[2]))
}

# exact segment-by-segment solution of dy/dt = -k(t) y for a step k
step_decay <- function(t, y0, bps, ks) {
  vapply(t, function(tt) {
    y <- y0
    for (i in seq_along(bps)) {
      a <- bps[i]
      b <- if (i < length(bps)) bps[i + 1L] else Inf
      if (tt <= a) break
      y <- y * exp(-ks[i] * (min(tt, b) - a))
    }
    y
  }, 0)
}

tv_fd <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}

# --- validation, no solver needed -----------------------------------

test_that("tv and tv_break are validated against each other", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3)
  expect_error(
    frm_ode(tv_decay, init = list(100), times = tt, parms = list(0.2),
            tv_break = c("a", "a", "b")),
    "`tv_break` was given but `tv` was not"
  )
  expect_error(
    frm_ode(tv_decay, init = list(100), times = tt, parms = list(0.2),
            tv = list(c(0.2, 0.2, 0.5)), tv_break = c("a", "b")),
    "it is one value per observation"
  )
  # a tv column that contradicts the declared blocks is a different
  # model from the one the breaks describe
  expect_error(
    frm_ode(tv_decay, init = list(100), times = tt,
            tv = list(c(0.2, 0.3, 0.5)), tv_break = c("a", "a", "b")),
    "changes inside a block"
  )
  # a step function of time cannot take two values at one time
  expect_error(
    frm_ode(tv_decay, init = list(100), times = c(1, 1, 3),
            tv = list(c(0.2, 0.5, 0.5))),
    "disagree about a `tv` value"
  )
  # no constants at all is a mistake, not an empty system
  expect_error(
    frm_ode(tv_decay, init = list(100), times = tt),
    "has no parameters at all"
  )
})

test_that("an estimated tv column needs its change points as data", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3)
  tp <- function() {
    RTMB::MakeTape(function(th) {
      "c" <- RTMB::ADoverload("c")
      sum(frm_ode(tv_decay, init = list(100), times = tt,
                  tv = list(exp(th[1] + th[2] * c(0, 0, 1)))))
    }, c(log(0.2), 0.5))
  }
  expect_error(tp(), "Pass `tv_break`")
  expect_error(
    frm_ode(tv_decay, init = list(100), times = tt, parms = list(0.2),
            tv = list(c(0.2, 0.2, 0.5)),
            tv_break = RTMB::advector(c(1, 1, 2))),
    "`tv_break` is an estimated quantity"
  )
})

# --- numerical correctness ------------------------------------------

test_that("a piecewise-constant rate matches the segment-by-segment solution", {
  skip_if_not_installed("RTMBode")
  # k = 0.2 up to t = 3, then 0.5, then 0.1 from t = 6
  tt <- c(0, 1, 2, 3, 4, 5, 6, 7, 9, 12)
  k <- ifelse(tt < 3, 0.2, ifelse(tt < 6, 0.5, 0.1))
  got <- frm_ode(tv_decay, init = list(100), times = tt, tv = list(k),
                 atol = 1e-12, rtol = 1e-12)
  ref <- step_decay(tt, 100, c(0, 3, 6), c(0.2, 0.5, 0.1))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("LOCF: a change point takes effect from its own row forward", {
  skip_if_not_installed("RTMBode")
  # The row at t = 3 carries 0.5, so 0.5 is in force on [3, .). The
  # state is continuous there - a covariate moves the derivative, not
  # the state - so the observation AT the change point still reads what
  # the pre-change rate produced, and the difference appears only after.
  tt <- c(2, 3, 3.5)
  k <- c(0.2, 0.5, 0.5)
  got <- frm_ode(tv_decay, init = list(100), times = tt, tv = list(k),
                 atol = 1e-12, rtol = 1e-12)
  expect_equal(got[2], 100 * exp(-0.2 * 3), tolerance = 1e-9)
  expect_equal(got[3], 100 * exp(-0.2 * 3) * exp(-0.5 * 0.5),
               tolerance = 1e-9)
  # the alternative reading, in which the new value starts at the NEXT
  # row, is a different number at t = 3.5 and is not what is shipped
  expect_false(isTRUE(all.equal(
    got[[3]], 100 * exp(-0.2 * 3.5), tolerance = 1e-6)))
  # before the first row the first block's value reaches back to t0
  early <- frm_ode(tv_decay, init = list(100), times = c(1, 4),
                   tv = list(c(0.2, 0.5)), atol = 1e-12, rtol = 1e-12)
  expect_equal(early[1], 100 * exp(-0.2 * 1), tolerance = 1e-9)
})

test_that("a constant tv column is the same model as a constant parm", {
  skip_if_not_installed("RTMBode")
  tt <- c(0.5, 2, 5)
  a <- frm_ode(tv_decay, init = list(100), times = tt, tv = list(0.3),
               atol = 1e-12, rtol = 1e-12)
  b <- frm_ode(tv_decay, init = list(100), times = tt, parms = list(0.3),
               atol = 1e-12, rtol = 1e-12)
  expect_equal(a, b, tolerance = 1e-10)
})

test_that("tv change points are unioned with the event times", {
  skip_if_not_installed("RTMBode")
  # rate changes at t = 4, a dose lands at t = 6: both split the solve,
  # and the answer is the hand-chained solution
  tt <- c(1, 4, 5, 6, 8, 10)
  k <- ifelse(tt < 4, 0.2, 0.6)
  got <- frm_ode(tv_decay, init = list(100), times = tt, tv = list(k),
                 events = data.frame(time = 6, value = 100),
                 atol = 1e-12, rtol = 1e-12)
  ref <- vapply(tt, function(u) {
    y <- step_decay(min(u, 6), 100, c(0, 4), c(0.2, 0.6))
    if (u <= 6) return(y)
    (y + 100) * exp(-0.6 * (u - 6))
  }, 0)
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("each group gets its own change points", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 3, 5, 7)
  d <- data.frame(id = factor(rep(c("a", "b"), each = 4)),
                  time = rep(tt, 2))
  # a switches at 3, b at 5
  d$k <- c(0.2, 0.6, 0.6, 0.6, 0.2, 0.2, 0.6, 0.6)
  got <- frm_ode(tv_decay, init = list(100), times = d$time,
                 group = d$id, tv = list(d$k),
                 atol = 1e-12, rtol = 1e-12)
  expect_equal(got[d$id == "a"],
               step_decay(tt, 100, c(0, 3), c(0.2, 0.6)), tolerance = 1e-8)
  expect_equal(got[d$id == "b"],
               step_decay(tt, 100, c(0, 5), c(0.2, 0.6)), tolerance = 1e-8)
})

test_that("a time-varying input reaches a two-state system by position", {
  skip_if_not_installed("RTMBode")
  # parms = (ka, V), tv = (ke), so the dynamics reads ke at p[3]
  tt <- c(1, 2, 4, 6, 8)
  ke <- ifelse(tt < 4, 0.2, 0.4)
  got <- frm_ode(tv_pk, init = list(100, 0), times = tt,
                 parms = list(1, 10), tv = list(ke),
                 states = c("depot", "central"), output = "central",
                 atol = 1e-12, rtol = 1e-12)
  # solve the same system by hand as two chained ordinary solves
  a <- frm_ode(tv_pk, init = list(100, 0), times = c(1, 2, 4),
               parms = list(1, 10), tv = list(0.2), output = 2L,
               atol = 1e-12, rtol = 1e-12)
  dep4 <- 100 * exp(-4)
  b <- frm_ode(tv_pk, init = list(dep4, a[3]), times = c(6, 8),
               parms = list(1, 10), tv = list(0.4), t0 = 4, output = 2L,
               atol = 1e-12, rtol = 1e-12)
  expect_equal(got, c(a, b), tolerance = 1e-8)
})

# --- the adjoint ----------------------------------------------------

test_that("an estimated time-varying rate is differentiated exactly", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  tt <- c(1, 2, 3, 5, 8)
  phase <- c("a", "a", "b", "b", "b")     # the change point is t = 3
  th0 <- c(log(0.2), 0.9)
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    k <- exp(th[1] + th[2] * (phase == "b"))
    sum(frm_ode(tv_decay, init = list(100), times = tt, tv = list(k),
                tv_break = phase, atol = 1e-11, rtol = 1e-11))
  }, th0)
  ref <- function(th) {
    sum(step_decay(tt, 100, c(0, 3), c(exp(th[1]), exp(th[1] + th[2]))))
  }
  expect_equal(tp(th0), ref(th0), tolerance = 1e-7)
  expect_equal(as.numeric(tp$jacfun()(th0)), tv_fd(ref, th0),
               tolerance = 1e-6)
})

# --- inside a model -------------------------------------------------

test_that("a tv nonlinear parameter is exempt from the constancy guard", {
  skip_if_not_installed("RTMBode")
  set.seed(21)
  d <- data.frame(id = factor(rep(1:3, each = 4)),
                  time = rep(c(1, 4, 8, 12), 3))
  d$phase <- factor(ifelse(d$time <= 4, "early", "late"))
  d$conc <- abs(stats::rnorm(nrow(d), 40, 5))
  st <- list(beta = c(log(0.2), 0))

  # in `parms` the same within-group covariate is refused by name
  refused <- bf(
    conc ~ frm_ode(tv_decay, init = list(100), times = time,
                   parms = list(exp(lk)), group = id, output = 1L),
    lk ~ 1 + phase, nl = TRUE)
  expect_error(
    frm(refused + gaussian(), data = d, dry_run = "frame", start = st),
    "not constant within 'id'"
  )
  # in `tv` it is the point of the feature
  allowed <- bf(
    conc ~ frm_ode(tv_decay, init = list(100), times = time,
                   tv = list(exp(lk)), tv_break = phase, group = id,
                   output = 1L),
    lk ~ 1 + phase, nl = TRUE)
  expect_s3_class(
    frm(allowed + gaussian(), data = d, dry_run = "frame", start = st),
    "frmtmb_frame")
})

test_that("a fit recovers a clearance that changes at a known time", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  set.seed(77)
  n_id <- 10
  tt <- c(0.5, 2, 4, 8, 11.9, 13, 16, 20)
  d <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                  time = rep(tt, n_id))
  d$phase <- factor(ifelse(d$time < 12, "early", "late"))
  b_i <- stats::rnorm(n_id, 0, 0.25)[as.integer(d$id)]
  k_early <- exp(log(0.15) + b_i)
  k_late <- exp(log(0.15) + log(2) + b_i)
  mu <- vapply(seq_len(nrow(d)), function(i) {
    step_decay(d$time[i], 100, c(0, 12), c(k_early[i], k_late[i]))
  }, 0)
  d$conc <- mu + stats::rnorm(nrow(d), 0, 0.6)

  fit <- frm(
    bf(conc ~ frm_ode(tv_decay, init = list(100), times = time,
                      tv = list(exp(lk)), tv_break = phase, group = id,
                      output = 1L),
       lk ~ 1 + phase + (1 | id), nl = TRUE) + gaussian(),
    data = d, start = list(beta = c(log(0.2), 0)))

  fx <- unlist(fixef(fit))
  expect_equal(unname(fx[["lk.(Intercept)"]]), log(0.15), tolerance = 0.15)
  expect_equal(unname(fx[["lk.phaselate"]]), log(2), tolerance = 0.3)
  # sigma comes back high (about 0.9 against a truth of 0.6): with ten
  # subjects and eight points the per-subject rate is itself estimated
  # with error, and the residual absorbs it. The bound is the assertion
  # that matters - the step model has not turned the noise into signal.
  expect_lt(unname(exp(fx[["sigma.(Intercept)"]])), 1.3)
  expect_true(all(is.finite(
    summary(fit)$coefficients$lk[, "Std. Error"])))

  # the shift is real: the same model without it fits worse
  flat <- frm(
    bf(conc ~ frm_ode(tv_decay, init = list(100), times = time,
                      parms = list(exp(lk)), group = id, output = 1L),
       lk ~ 1 + (1 | id), nl = TRUE) + gaussian(),
    data = d, start = list(beta = log(0.2)))
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(flat)) + 10)

  # predict() re-solves the step function on newdata
  nd <- d[d$id %in% c("1", "2"), ]
  expect_equal(predict(fit, newdata = nd),
               unname(predict(fit)[d$id %in% c("1", "2")]),
               tolerance = 1e-6, ignore_attr = TRUE)
})
