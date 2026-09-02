# frm_ode(events = ): repeated doses and infusions.
#
# The reference throughout is the analytic multi-dose superposition for
# the one-compartment oral model, which is exact, so a mismatch is the
# event machinery and nothing else. dev/ode/probeH1..H4 establish what
# the backend does and does not support; the two facts that shape this
# file are:
#
#   - deSolve's own `events` argument is unusable here. It errors on the
#     automatic-differentiation path (RTMBode hands deSolve an unnamed
#     state vector), and if that were fixed it would give a WRONG
#     gradient for "replace" and "multiply", because the event jumps the
#     state without jumping the sensitivity block RTMBode integrates
#     alongside it. frm_ode() splits the solve at the event times
#     instead, which is exact for every method.
#   - The dose amount may therefore depend on a parameter, through
#     `event_scale`.
#
# As in test-ode.R, every test that fails a solve on purpose lives at
# the end, and its warnings are captured rather than let into the test
# record.

pk_dyn2 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

decay1 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}

pk_analytic2 <- function(t, ka, ke, V, D) {
  D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
}

# C(t) for a sequence of oral bolus doses, by superposition
multi_dose <- function(t, ka, ke, V, amt, dose_times) {
  vapply(seq_along(t), function(i) {
    keep <- dose_times <= t[i]
    if (!any(keep)) return(0)
    u <- t[i] - dose_times[keep]
    sum(amt[keep] * ka / (V * (ka - ke)) *
          (exp(-ke * u) - exp(-ka * u)))
  }, 0)
}

# A(t) for a constant-rate infusion of `amt` over `dur` into a
# one-compartment system with elimination `ke`
infusion1 <- function(t, amt, dur, ke) {
  R <- amt / dur
  ifelse(t <= dur, R / ke * (1 - exp(-ke * t)),
         R / ke * (1 - exp(-ke * dur)) * exp(-ke * (t - dur)))
}

central_fd <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}

# --- validation, no solver needed -----------------------------------

test_that("the events table is validated", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(rep(1:2, each = 3)),
                  time = rep(c(1, 2, 3), 2))
  call_ev <- function(ev, ...) {
    frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id,
            states = c("depot", "central"), output = "central",
            events = ev, ...)
  }
  expect_error(call_ev(list(time = 1, value = 1)), "must be a data.frame")
  expect_error(call_ev(data.frame()[0, ]), "no rows")
  expect_error(call_ev(data.frame(time = 1)), "missing the value column")
  expect_error(call_ev(data.frame(value = 1)), "missing the time column")
  expect_error(call_ev(data.frame(time = 1, value = 1, evid = 1)),
               "unknown column: evid")
  expect_error(call_ev(data.frame(time = 1, value = 1, evid = 1)),
               "does not read NONMEM records")
  expect_error(call_ev(data.frame(time = NA_real_, value = 1)),
               "must be finite and numeric")
  expect_error(call_ev(data.frame(time = 1, value = NA_real_)),
               "must be finite and numeric")
  # two states, so a compartment has to be named
  expect_error(call_ev(data.frame(time = 1, value = 1)),
               "no `state` column")
  expect_error(call_ev(data.frame(time = 1, value = 1, state = "gut")),
               "not in `states`")
  expect_error(call_ev(data.frame(time = 1, value = 1, state = 7L)),
               "index states 1 to 2")
  expect_error(
    call_ev(data.frame(time = 1, value = 1, state = 1L,
                       method = "bolus")),
    "unknown method"
  )
  expect_error(
    call_ev(data.frame(time = 1, value = 1, state = 1L, duration = -1)),
    "not be negative|not negative"
  )
  expect_error(
    call_ev(data.frame(time = 1, value = 1, state = 1L, duration = 2,
                       method = "replace")),
    "is positive on a row whose method is not"
  )
  expect_error(
    call_ev(data.frame(time = 1, value = 1, state = 1L, group = "9")),
    "not in `group`"
  )
  # an event before t0 is a data error, not something to integrate over
  expect_error(
    call_ev(data.frame(time = 0.5, value = 1, state = 1L), t0 = 1),
    "before t0"
  )
  # two rows on one state at one instant compose only as additions
  expect_error(
    call_ev(data.frame(time = c(1, 1), value = c(1, 2), state = 1L,
                       method = c("add", "replace"))),
    "ambiguous"
  )
  expect_silent(call_ev(data.frame(time = c(1, 1), value = c(1, 2),
                                   state = 1L)))
})

test_that("event_scale is refused where scaling would change a meaning", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3)
  expect_error(
    frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
            event_scale = 0.5),
    "`event_scale` was given but `events` was not"
  )
  expect_error(
    frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
            events = data.frame(time = 1.5, value = 50,
                                method = "replace"),
            event_scale = 0.5),
    "applies only to"
  )
})

test_that("a function-valued events table is called", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3)
  sched <- function() data.frame(time = 1.5, value = 100)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(0.2),
                 events = sched, atol = 1e-10, rtol = 1e-10)
  expect_equal(got, c(0, 100 * exp(-0.2 * 0.5), 100 * exp(-0.2 * 1.5)),
               tolerance = 1e-7)
  expect_error(
    frm_ode(decay1, init = list(0), times = tt, parms = list(0.2),
            events = function() 1),
    "must return a data.frame"
  )
})

test_that("a bare events table in a formula is the inline table", {
  skip_if_not_installed("RTMBode")
  # dev/ode-feasibility.md section 9.6: `events = doses` used to reach
  # model.frame() as a request for a column called `doses` and die with
  # "invalid type (list)". It is an argument of frm_ode(), not a column,
  # so it now resolves in the formula environment - and it must give
  # back exactly the model the inline data.frame gives.
  set.seed(4)
  n_id <- 4
  tt <- c(0.5, 2, 6, 11.9, 14, 20, 26)
  d <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                  time = rep(tt, n_id))
  d$conc <- abs(stats::rnorm(nrow(d), 5, 1))
  doses <- data.frame(time = c(12, 24), state = "depot", value = 100)

  form_sym <- bf(
    conc ~ frm_ode(pk_dyn2, init = list(100, 0), times = time,
                   parms = list(exp(lka), exp(lke), exp(lV)), group = id,
                   states = c("depot", "central"), output = "central",
                   events = doses),
    lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE)
  form_inl <- bf(
    conc ~ frm_ode(pk_dyn2, init = list(100, 0), times = time,
                   parms = list(exp(lka), exp(lke), exp(lV)), group = id,
                   states = c("depot", "central"), output = "central",
                   events = data.frame(time = c(12, 24), state = "depot",
                                       value = 100)),
    lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE)
  st <- list(beta = c(0, log(0.25), log(8)))

  fr <- frm(form_sym + gaussian(), data = d, dry_run = "frame",
            start = st)
  expect_false("doses" %in% names(fr$linpreds[["conc.mu"]]$data_list))
  expect_setequal(names(fr$linpreds[["conc.mu"]]$data_list),
                  c("time", "id"))
  expect_true("doses" %in% fr$linpreds[["conc.mu"]]$nl_lexical)

  f_sym <- frm(form_sym + gaussian(), data = d, start = st)
  f_inl <- frm(form_inl + gaussian(), data = d, start = st)
  expect_equal(as.numeric(logLik(f_sym)), as.numeric(logLik(f_inl)),
               tolerance = 1e-12)
  expect_equal(unlist(fixef(f_sym)), unlist(fixef(f_inl)),
               tolerance = 1e-12)
})

# --- numerical correctness ------------------------------------------

test_that("repeated doses match the analytic superposition", {
  skip_if_not_installed("RTMBode")
  dose_t <- c(0, 12, 24, 36)
  obs <- c(0.5, 2, 6, 11, 13, 18, 25, 30, 37, 42, 48)
  d <- data.frame(id = factor(rep(c("a", "b"), each = length(obs))),
                  time = rep(obs, 2))
  ref <- multi_dose(d$time, 1, 0.2, 10, rep(100, 4), dose_t)

  # first dose as the initial condition, the rest as events
  got <- frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id,
                 states = c("depot", "central"), output = "central",
                 events = data.frame(time = dose_t[-1], state = "depot",
                                     value = 100),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got, ref, tolerance = 1e-8)

  # every dose as an event, starting from an empty depot: the same curve
  got2 <- frm_ode(pk_dyn2, init = list(0, 0), times = d$time,
                  parms = list(1, 0.2, 10), group = d$id,
                  states = c("depot", "central"), output = "central",
                  events = data.frame(time = dose_t, state = "depot",
                                      value = 100),
                  atol = 1e-10, rtol = 1e-10)
  expect_equal(got2, ref, tolerance = 1e-8)
})

test_that("a per-group schedule doses only its own group", {
  skip_if_not_installed("RTMBode")
  obs <- c(1, 6, 13, 20, 30)
  d <- data.frame(id = factor(rep(c("a", "b"), each = length(obs))),
                  time = rep(obs, 2))
  ev <- data.frame(group = c("a", "a", "b"), time = c(12, 24, 12),
                   state = "depot", value = 100)
  got <- frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id,
                 states = c("depot", "central"), output = "central",
                 events = ev, atol = 1e-10, rtol = 1e-10)
  expect_equal(got[d$id == "a"],
               multi_dose(obs, 1, 0.2, 10, rep(100, 3), c(0, 12, 24)),
               tolerance = 1e-8)
  expect_equal(got[d$id == "b"],
               multi_dose(obs, 1, 0.2, 10, rep(100, 2), c(0, 12)),
               tolerance = 1e-8)
})

test_that("a group with no events takes the plain single-solve path", {
  skip_if_not_installed("RTMBode")
  obs <- c(1, 6, 13, 20)
  d <- data.frame(id = factor(rep(c("a", "b"), each = length(obs))),
                  time = rep(obs, 2))
  got <- frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id,
                 states = c("depot", "central"), output = "central",
                 events = data.frame(group = "a", time = 12,
                                     state = "depot", value = 100),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got[d$id == "b"],
               pk_analytic2(obs, 1, 0.2, 10, 100), tolerance = 1e-8)
})

test_that("an observation at a dose time reads the pre-dose value", {
  skip_if_not_installed("RTMBode")
  # the trough, matching deSolve's own reading at an event time
  tt <- c(0, 2, 4)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 events = data.frame(time = c(0, 2), value = 100),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got[1], 100)                       # init, before the t0 dose
  expect_equal(got[2], 200 * exp(-0.2 * 2), tolerance = 1e-8)
  expect_equal(got[3], (200 * exp(-0.2 * 2) + 100) * exp(-0.2 * 2),
               tolerance = 1e-8)
})

test_that("duplicate times, ragged groups and unsorted rows still work", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(c(1, 1, 1, 1, 2, 2, 2)),
                  time = c(4, 1, 1, 0, 8, 0.5, 2))
  got <- frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
                 parms = list(1, 0.2, 10), group = d$id, output = 2L,
                 events = data.frame(time = 2, value = 100, state = 1L),
                 atol = 1e-10, rtol = 1e-10)
  ref <- multi_dose(d$time, 1, 0.2, 10, c(100, 100), c(0, 2))
  expect_equal(got, ref, tolerance = 1e-8)
  expect_equal(got[2], got[3])                    # the duplicate time
  expect_equal(got[4], 0)                         # C(0) = 0
})

test_that("an unsorted events table is put in order", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 5, 9, 15)
  a <- frm_ode(decay1, init = list(0), times = tt, parms = list(0.2),
               events = data.frame(time = c(8, 0, 4), value = 100),
               atol = 1e-10, rtol = 1e-10)
  b <- frm_ode(decay1, init = list(0), times = tt, parms = list(0.2),
               events = data.frame(time = c(0, 4, 8), value = 100),
               atol = 1e-10, rtol = 1e-10)
  expect_equal(a, b)
})

test_that("t0 shifts the origin of a dosing solve", {
  skip_if_not_installed("RTMBode")
  tt <- c(2.5, 5, 9)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 t0 = 1,
                 events = data.frame(time = 4, value = 100),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got,
               c(100 * exp(-0.2 * 1.5),
                 (100 * exp(-0.2 * 3) + 100) * exp(-0.2 * 1),
                 (100 * exp(-0.2 * 3) + 100) * exp(-0.2 * 5)),
               tolerance = 1e-8)
})

test_that("replace and multiply are exact", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3, 4)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.25),
                 events = data.frame(time = 2, value = 50,
                                     method = "replace"),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got, c(100 * exp(-0.25), 100 * exp(-0.5),
                      50 * exp(-0.25), 50 * exp(-0.5)),
               tolerance = 1e-8)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.25),
                 events = data.frame(time = 2, value = 2,
                                     method = "multiply"),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got, c(100 * exp(-0.25), 100 * exp(-0.5),
                      200 * exp(-0.75), 200 * exp(-1)),
               tolerance = 1e-8)
})

test_that("events beyond the last observation change nothing", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3)
  a <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
               events = data.frame(time = 10, value = 500))
  expect_equal(a, 100 * exp(-0.2 * tt), tolerance = 1e-7)
})

test_that("a group observed only at t0 reads init, dose or no dose", {
  skip_if_not_installed("RTMBode")
  # no segment to integrate, and the t0 dose is still after the reading
  got <- frm_ode(decay1, init = list(5), times = c(0, 0), parms = list(0.2),
                 events = data.frame(time = 0, value = 100))
  expect_equal(got, c(5, 5))
})

test_that("an infusion still running at the last observation works", {
  skip_if_not_installed("RTMBode")
  # the window is truncated at the last observation, and the rate is
  # active over every segment up to it
  tt <- c(1, 2)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(0.2),
                 events = data.frame(time = 0, value = 80, duration = 10),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got, (80 / 10) / 0.2 * (1 - exp(-0.2 * tt)),
               tolerance = 1e-8)
})

test_that("an infusion matches the analytic constant-rate solution", {
  skip_if_not_installed("RTMBode")
  tt <- c(0, 1, 2, 4, 4.5, 6, 10, 20)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(0.25),
                 events = data.frame(time = 0, value = 80, duration = 4),
                 atol = 1e-10, rtol = 1e-10)
  expect_equal(got, infusion1(tt, 80, 4, 0.25), tolerance = 1e-7)

  # repeated infusions superpose
  tt2 <- c(1, 3, 5, 9, 13, 15, 20)
  got2 <- frm_ode(decay1, init = list(0), times = tt2, parms = list(0.25),
                  events = data.frame(time = c(0, 12), value = 80,
                                      duration = 4),
                  atol = 1e-10, rtol = 1e-10)
  expect_equal(got2,
               infusion1(tt2, 80, 4, 0.25) +
                 ifelse(tt2 > 12, infusion1(tt2 - 12, 80, 4, 0.25), 0),
               tolerance = 1e-7)
})

test_that("a bolus and an infusion into different states coexist", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 3, 5)
  got <- frm_ode(pk_dyn2, init = list(0, 0), times = tt,
                 parms = list(1, 0.2, 10),
                 states = c("depot", "central"),
                 events = data.frame(
                   time = c(0, 0), state = c("depot", "central"),
                   value = c(100, 40), duration = c(0, 4)),
                 atol = 1e-10, rtol = 1e-10)
  expect_identical(dim(got), c(3L, 2L))
  # the depot only ever sees the bolus
  expect_equal(got[, 1], 100 * exp(-tt), tolerance = 1e-7)
  # the central compartment is the oral curve plus the infusion
  expect_equal(got[, 2],
               pk_analytic2(tt, 1, 0.2, 10, 100) +
                 infusion1(tt, 40, 4, 0.2),
               tolerance = 1e-7)
})

# --- the adjoint ----------------------------------------------------

test_that("the gradient through repeated doses matches finite differences", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  dose_t <- c(0, 12, 24)
  obs <- c(1, 6, 13, 20, 26, 34)
  th0 <- c(0, log(0.2), log(10))
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    n <- length(obs)
    sum(frm_ode(pk_dyn2, init = list(0, 0), times = obs,
                parms = list(rep(exp(th[1]), n), rep(exp(th[2]), n),
                             rep(exp(th[3]), n)),
                states = c("depot", "central"), output = "central",
                events = data.frame(time = dose_t, state = "depot",
                                    value = 100),
                atol = 1e-10, rtol = 1e-10))
  }, th0)
  ref <- function(th)
    sum(multi_dose(obs, exp(th[1]), exp(th[2]), exp(th[3]),
                   rep(100, 3), dose_t))
  expect_equal(as.numeric(tp$jacfun()(th0)), central_fd(ref, th0),
               tolerance = 1e-6)
})

test_that("the gradient through replace and multiply is exact", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # this is the case a deSolve events table gets wrong: the state jumps
  # and the sensitivity block RTMBode carries alongside it does not
  tt <- c(1, 2, 3, 4)
  for (m in c("replace", "multiply")) {
    val <- if (m == "replace") 50 else 2
    tp <- RTMB::MakeTape(function(th) {
      "c" <- RTMB::ADoverload("c")
      sum(frm_ode(decay1, init = list(100), times = tt,
                  parms = list(exp(th[1])),
                  events = data.frame(time = 2, value = val, method = m),
                  atol = 1e-10, rtol = 1e-10))
    }, log(0.25))
    ref <- function(th) {
      k <- exp(th[1])
      if (m == "replace") {
        sum(c(100 * exp(-k), 100 * exp(-2 * k),
              50 * exp(-k), 50 * exp(-2 * k)))
      } else {
        sum(c(100 * exp(-k), 100 * exp(-2 * k),
              200 * exp(-3 * k), 200 * exp(-4 * k)))
      }
    }
    expect_equal(as.numeric(tp$jacfun()(log(0.25))),
                 central_fd(ref, log(0.25)), tolerance = 1e-6,
                 info = m)
  }
})

test_that("event_scale is an estimated dose multiplier", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  dose_t <- c(0, 12, 24)
  obs <- c(1, 6, 13, 20, 26, 34)
  th0 <- c(0, log(0.2), log(10), 0.4)
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    n <- length(obs)
    sum(frm_ode(pk_dyn2, init = list(0, 0), times = obs,
                parms = list(rep(exp(th[1]), n), rep(exp(th[2]), n),
                             rep(exp(th[3]), n)),
                states = c("depot", "central"), output = "central",
                events = data.frame(time = dose_t, state = "depot",
                                    value = 100),
                event_scale = rep(1 / (1 + exp(-th[4])), n),
                atol = 1e-10, rtol = 1e-10))
  }, th0)
  ref <- function(th)
    sum(multi_dose(obs, exp(th[1]), exp(th[2]), exp(th[3]),
                   rep(1 / (1 + exp(-th[4])) * 100, 3), dose_t))
  expect_equal(tp(th0), ref(th0), tolerance = 1e-7)
  expect_equal(as.numeric(tp$jacfun()(th0)), central_fd(ref, th0),
               tolerance = 1e-6)
})

test_that("an estimated infusion rate is differentiated exactly", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  tt <- c(1, 3, 5, 9, 13, 15, 20)
  x0 <- c(log(0.25), 0.3)
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    sum(frm_ode(decay1, init = list(0), times = tt,
                parms = list(exp(th[1])),
                events = data.frame(time = c(0, 12), value = 80,
                                    duration = 4),
                event_scale = 1 / (1 + exp(-th[2])),
                atol = 1e-10, rtol = 1e-10))
  }, x0)
  ref <- function(th) {
    amt <- 1 / (1 + exp(-th[2])) * 80
    k <- exp(th[1])
    sum(infusion1(tt, amt, 4, k) +
          ifelse(tt > 12, infusion1(tt - 12, amt, 4, k), 0))
  }
  expect_equal(tp(x0), ref(x0), tolerance = 1e-7)
  expect_equal(as.numeric(tp$jacfun()(x0)), central_fd(ref, x0),
               tolerance = 1e-6)
})

# --- inside a model -------------------------------------------------

test_that("a repeated-dosing population fit recovers the truth", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  set.seed(11)
  n_id <- 12
  dose_t <- c(0, 12, 24)
  tt <- c(0.5, 1, 2, 4, 8, 11.9, 13, 16, 23.9, 26, 30, 36)
  d <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                  time = rep(tt, n_id))
  ka <- exp(stats::rnorm(n_id, 0, 0.30))[as.integer(d$id)]
  ke <- exp(stats::rnorm(n_id, log(0.2), 0.25))[as.integer(d$id)]
  mu <- vapply(seq_len(nrow(d)), function(i)
    multi_dose(d$time[i], ka[i], ke[i], 10, rep(100, 3), dose_t), 0)
  d$conc <- mu + stats::rnorm(nrow(d), 0, 0.30)

  fit <- frm(
    bf(conc ~ frm_ode(pk_dyn2, init = list(100, 0), times = time,
                      parms = list(exp(lka), exp(lke), exp(lV)),
                      group = id, states = c("depot", "central"),
                      output = "central",
                      events = data.frame(time = c(12, 24),
                                          state = "depot", value = 100)),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = d, start = list(beta = c(0, log(0.25), log(8))))

  expect_s3_class(fit, "frmtmb_fit")
  fx <- unlist(fixef(fit))
  expect_equal(unname(fx[["lka.(Intercept)"]]), 0, tolerance = 0.4)
  expect_equal(unname(fx[["lke.(Intercept)"]]), log(0.2), tolerance = 0.3)
  expect_equal(unname(fx[["lV.(Intercept)"]]), log(10), tolerance = 0.3)
  expect_equal(unname(exp(fx[["sigma.(Intercept)"]])), 0.30,
               tolerance = 0.1)

  # the fitted curve is the multi-dose curve, not the single-dose one
  pr <- stats::predict(fit)
  single <- pk_analytic2(d$time, exp(fx[["lka.(Intercept)"]]),
                        exp(fx[["lke.(Intercept)"]]),
                        exp(fx[["lV.(Intercept)"]]), 100)
  expect_gt(max(abs(pr - single)), 1)
  expect_lt(stats::sd(pr - d$conc), 0.5)
})

test_that("a dosing model rejects a within-group dynamics covariate", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(rep(1:3, each = 4)),
                  time = rep(c(1, 6, 13, 20), 3),
                  x = rnorm(12), conc = rnorm(12, 5))
  form <- bf(
    conc ~ frm_ode(pk_dyn2, init = list(100, 0), times = time,
                   parms = list(exp(lka), exp(lke), exp(lV)),
                   group = id, output = 2L,
                   events = data.frame(time = 12, value = 100,
                                       state = 1L)),
    lka ~ 1 + x, lke ~ 1, lV ~ 1, nl = TRUE)
  expect_error(
    frm(form + gaussian(), data = d,
        start = list(beta = c(0, 0, log(0.25), log(8)))),
    "not constant within"
  )
})

test_that("event_scale must be constant within a group", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(rep(1:2, each = 3)),
                  time = rep(c(1, 6, 13), 2))
  expect_error(
    frm_ode(pk_dyn2, init = list(100, 0), times = d$time,
            parms = list(1, 0.2, 10), group = d$id, output = 2L,
            events = data.frame(time = 12, value = 100, state = 1L),
            event_scale = d$time),
    "`event_scale` column 1 is not constant within group"
  )
})

# --- failures, last: a bad solve can poison later tapes -------------

test_that("a failed dosing solve becomes a penalty naming the group", {
  skip_if_not_installed("RTMBode")
  # y' = y^2 runs away in finite time. At y(0) = 0.01 it is placid over
  # [0, 3]; the dose at t = 2 makes the NEXT segment blow up, so the
  # failure is specifically in the segmented path.
  blow <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(p[1] * y[1] * y[1]))
  }
  d <- data.frame(id = factor(rep(c("a", "b"), each = 2)),
                  time = rep(c(1, 3), 2))
  ev <- data.frame(group = "b", time = 2, value = 1e6, state = 1L)
  w <- capture_warnings(
    got <- frm_ode(blow, init = list(0.01), times = d$time,
                   parms = list(c(0, 0, 1, 1)), group = d$id,
                   events = ev, penalty = 1234)
  )
  expect_true(any(grepl("the solve failed for 1 of 2 groups", w)))
  expect_true(any(grepl("\\(b\\)", w)))
  expect_equal(unname(got[d$id == "b"]), c(1234, 1234))
  expect_equal(unname(got[d$id == "a"]), c(0.01, 0.01))
  fl <- frm_ode_failures()
  expect_identical(fl$groups, "b")
  expect_identical(fl$penalty, 1234)

  expect_error(
    suppressWarnings(
      frm_ode(blow, init = list(0.01), times = d$time,
              parms = list(c(0, 0, 1, 1)), group = d$id,
              events = ev, on_error = "error")),
    "failed to solve group 'b'"
  )
})
