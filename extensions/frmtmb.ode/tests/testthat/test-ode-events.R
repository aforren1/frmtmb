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

# --- ii / addl: compact repetition ----------------------------------

test_that("ii and addl expand to the hand-written rows exactly", {
  skip_if_not_installed("RTMBode")
  obs <- c(1, 6, 13, 20, 26, 34, 40, 50)
  compact <- frm_ode(decay1, init = list(0), times = obs,
                     parms = list(0.2),
                     events = data.frame(time = 0, value = 100, ii = 12,
                                         addl = 3),
                     atol = 1e-12, rtol = 1e-12)
  written <- frm_ode(decay1, init = list(0), times = obs,
                     parms = list(0.2),
                     events = data.frame(time = c(0, 12, 24, 36),
                                         value = 100),
                     atol = 1e-12, rtol = 1e-12)
  expect_identical(compact, written)

  # the expansion keeps every other column of the row it came from
  ev <- data.frame(group = "a", time = 2, state = 1L, value = 60,
                   method = "add", duration = 3, ii = 10, addl = 2)
  a <- frm_ode(decay1, init = list(0), times = obs, parms = list(0.2),
               group = rep("a", length(obs)), events = ev,
               atol = 1e-12, rtol = 1e-12)
  b <- frm_ode(decay1, init = list(0), times = obs, parms = list(0.2),
               events = data.frame(time = c(2, 12, 22), value = 60,
                                   duration = 3),
               atol = 1e-12, rtol = 1e-12)
  expect_identical(a, b)
})

test_that("ii and addl are validated", {
  skip_if_not_installed("RTMBode")
  call_ev <- function(ev) {
    frm_ode(decay1, init = list(0), times = c(1, 2), parms = list(0.2),
            events = ev)
  }
  expect_error(call_ev(data.frame(time = 0, value = 1, addl = 2)),
               "is positive on a row whose `ii` is not")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = -1)),
               "`events\\$ii` must be finite and not negative")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = 5,
                                  addl = 1.5)),
               "must be a whole number")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = "x")),
               "`events\\$ii` must be numeric")
})

# --- reset events (NONMEM EVID 3, rxode2 evid 3) --------------------

test_that("a reset zeroes every compartment", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3, 4)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 events = data.frame(time = 2, value = 0,
                                     method = "reset"),
                 atol = 1e-12, rtol = 1e-12)
  # the observation at the reset time is still the pre-event reading
  expect_equal(got, c(100 * exp(-0.2), 100 * exp(-0.4), 0, 0),
               tolerance = 1e-9)

  # every state, not the one a `state` column would name
  m <- frm_ode(pk_dyn2, init = list(100, 0), times = tt,
               parms = list(1, 0.2, 10),
               states = c("depot", "central"),
               events = data.frame(time = 2, value = 0,
                                   method = "reset"),
               atol = 1e-12, rtol = 1e-12)
  expect_equal(unname(m[3:4, 1]), c(0, 0))
  expect_equal(unname(m[3:4, 2]), c(0, 0))

  # a reset needs no compartment even when the system has several
  expect_silent(frm_ode(pk_dyn2, init = list(100, 0), times = tt,
                        parms = list(1, 0.2, 10),
                        events = data.frame(time = 2, value = 0,
                                            method = "reset")))
})

test_that("a reset beside a dose is EVID 4, in that order", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 2, 3, 4)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 events = data.frame(time = c(2, 2), value = c(0, 50),
                                     method = c("reset", "add")),
                 atol = 1e-12, rtol = 1e-12)
  expect_equal(got, c(100 * exp(-0.2), 100 * exp(-0.4),
                      50 * exp(-0.2), 50 * exp(-0.4)),
               tolerance = 1e-9)
  # the row order in the table does not decide it
  rev <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 events = data.frame(time = c(2, 2), value = c(50, 0),
                                     method = c("add", "reset")),
                 atol = 1e-12, rtol = 1e-12)
  expect_identical(got, rev)
})

test_that("a reset to a non-zero level sets every state to it", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 3)
  got <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                 events = data.frame(time = 2, value = 7,
                                     method = "reset"),
                 atol = 1e-12, rtol = 1e-12)
  expect_equal(got[2], 7 * exp(-0.2), tolerance = 1e-9)
})

test_that("event_scale does not scale a reset", {
  skip_if_not_installed("RTMBode")
  tt <- c(1, 3, 5)
  ev <- data.frame(time = c(2, 4), value = c(7, 100),
                   method = c("reset", "add"))
  half <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                  events = ev, event_scale = 0.5,
                  atol = 1e-12, rtol = 1e-12)
  full <- frm_ode(decay1, init = list(100), times = tt, parms = list(0.2),
                  events = ev, atol = 1e-12, rtol = 1e-12)
  # the reset level is the same either way, the dose is halved
  expect_equal(half[2], full[2], tolerance = 1e-10)
  expect_equal(half[3], (7 * exp(-0.4) + 50) * exp(-0.2),
               tolerance = 1e-9)
  expect_equal(full[3], (7 * exp(-0.4) + 100) * exp(-0.2),
               tolerance = 1e-9)
})

test_that("the gradient through a reset is exact", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  tt <- c(1, 2, 3, 4)
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    sum(frm_ode(decay1, init = list(100), times = tt,
                parms = list(exp(th[1])),
                events = data.frame(time = c(2, 2), value = c(0, 50),
                                    method = c("reset", "add")),
                atol = 1e-11, rtol = 1e-11))
  }, log(0.25))
  ref <- function(th) {
    k <- exp(th[1])
    sum(c(100 * exp(-k), 100 * exp(-2 * k), 50 * exp(-k), 50 * exp(-2 * k)))
  }
  expect_equal(tp(log(0.25)), ref(log(0.25)), tolerance = 1e-7)
  expect_equal(as.numeric(tp$jacfun()(log(0.25))),
               central_fd(ref, log(0.25)), tolerance = 1e-6)
})

# --- steady-state dosing --------------------------------------------

# The exact steady state of a one-compartment system under a bolus D
# every tau: the trough is a geometric series, D e^-k tau / (1 - e^-k tau).
ss_trough <- function(D, k, tau) D * exp(-k * tau) / (1 - exp(-k * tau))

test_that("a steady-state row matches the analytic superposition", {
  skip_if_not_installed("RTMBode")
  k <- 0.2
  tt <- c(0, 1, 3, 6, 11.999)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(k),
                 events = data.frame(time = 0, value = 100, ii = 12,
                                     ss = TRUE),
                 atol = 1e-12, rtol = 1e-12)
  tr <- ss_trough(100, k, 12)
  ref <- (tr + 100) * exp(-k * tt)
  ref[1] <- tr            # the record's own time reads the trough
  expect_equal(got, ref, tolerance = 1e-8)

  # the run-in replaces `init`, whatever it held
  got2 <- frm_ode(decay1, init = list(500), times = tt, parms = list(k),
                  events = data.frame(time = 0, value = 100, ii = 12,
                                      ss = TRUE),
                  atol = 1e-12, rtol = 1e-12)
  expect_equal(got2, ref, tolerance = 1e-8)
})

test_that("a steady-state row continues into ordinary doses", {
  skip_if_not_installed("RTMBode")
  k <- 0.2
  tt <- c(0, 6, 12, 18, 24)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(k),
                 events = data.frame(time = 0, value = 100, ii = 12,
                                     addl = 2, ss = TRUE),
                 atol = 1e-12, rtol = 1e-12)
  tr <- ss_trough(100, k, 12)
  # at steady state the cycle repeats, so every trough and every
  # mid-interval point is the same number
  expect_equal(got, c(tr, (tr + 100) * exp(-6 * k), tr,
                      (tr + 100) * exp(-6 * k), tr), tolerance = 1e-8)
})

test_that("a steady-state infusion matches the infinite superposition", {
  skip_if_not_installed("RTMBode")
  k <- 0.2
  dur <- 3
  tau <- 12
  tt <- c(0, 1, 3, 6, 11)
  got <- frm_ode(decay1, init = list(0), times = tt, parms = list(k),
                 events = data.frame(time = 0, value = 100,
                                     duration = dur, ii = tau, ss = TRUE),
                 n_ss = 30, atol = 1e-12, rtol = 1e-12)
  one <- function(t) {
    R <- 100 / dur
    ifelse(t <= dur, R / k * (1 - exp(-k * t)),
           R / k * (1 - exp(-k * dur)) * exp(-k * (t - dur)))
  }
  ref <- vapply(tt, function(t) sum(vapply(0:400, function(j)
    one(t + j * tau), 0)), 0)
  expect_equal(got, ref, tolerance = 1e-7)
})

test_that("n_ss controls the run-in and its error falls off geometrically", {
  skip_if_not_installed("RTMBode")
  k <- 0.2
  tr <- ss_trough(100, k, 12)
  # ss_extrapolate = FALSE: this measures the TRUNCATION law, and the
  # default sums the tail that law describes rather than dropping it
  # ss_tol = Inf silences the shortfall warning that is the point of
  # the next test; here the shortfall itself is what is measured
  err <- vapply(c(2L, 4L, 8L), function(n) {
    g <- frm_ode(decay1, init = list(0), times = 0, parms = list(k),
                 events = data.frame(time = 0, value = 100, ii = 12,
                                     ss = TRUE),
                 n_ss = n, ss_tol = Inf, ss_extrapolate = FALSE,
                 atol = 1e-12, rtol = 1e-12)
    abs(g - tr) / tr
  }, 0)
  # exactly exp(-n k tau) for linear kinetics
  expect_equal(err, exp(-c(2, 4, 8) * k * 12), tolerance = 1e-6)
  expect_error(
    frm_ode(decay1, init = list(0), times = 0, parms = list(k),
            events = data.frame(time = 0, value = 100, ii = 12,
                                ss = TRUE), n_ss = 0),
    "`n_ss` must be one whole number"
  )
})

test_that("an unconverged run-in warns on the numeric path", {
  skip_if_not_installed("RTMBode")
  # a drug whose half-life is long against its interval: three truncated
  # cycles are nowhere near the steady state, and only the numeric path
  # can say so
  expect_warning(
    frm_ode(decay1, init = list(0), times = c(0, 5), parms = list(0.005),
            events = data.frame(time = 0, value = 100, ii = 6,
                                ss = TRUE), n_ss = 3,
            ss_extrapolate = FALSE),
    "steady-state run-in cycles"
  )
  expect_silent(
    frm_ode(decay1, init = list(0), times = c(0, 5), parms = list(0.005),
            events = data.frame(time = 0, value = 100, ii = 6,
                                ss = TRUE), n_ss = 3, ss_tol = Inf,
            ss_extrapolate = FALSE)
  )
})

# --- the steady-state tail ------------------------------------------

test_that("the geometric tail is summed rather than dropped", {
  skip_if_not_installed("RTMBode")
  # k * ii = 0.12, a terminal half-life of about six dosing intervals,
  # which is where the truncated run-in is worst
  k <- 0.01
  tau <- 12
  tr <- ss_trough(100, k, tau)
  ev <- data.frame(time = 0, value = 100, ii = tau, ss = TRUE)
  got <- function(ext) {
    frm_ode(decay1, init = list(0), times = 0, parms = list(k),
            events = ev, n_ss = 20L, ss_tol = Inf, ss_extrapolate = ext,
            atol = 1e-10, rtol = 1e-10)
  }
  e_trunc <- abs(got(FALSE) - tr) / tr
  e_ext <- abs(got(TRUE) - tr) / tr
  # the shortfall the truncation leaves is exp(-n k tau), so this pins
  # the size of the defect and not only its direction
  expect_equal(as.numeric(e_trunc), exp(-20 * k * tau),
               tolerance = 1e-5)
  # the bar is a ratio to what the same run measured, never a constant
  expect_lt(as.numeric(e_ext), as.numeric(e_trunc) / 1e4)
})

test_that("a state with no steady state does not disable the tail", {
  skip_if_not_installed("RTMBode")
  # the second state is an AUC compartment: it grows without bound, so
  # its cycle-to-cycle change never shrinks. Read one shared ratio over
  # all states and it is that state that sets it, which stands the
  # correction down for the state that DOES settle.
  k <- 0.01
  tau <- 12
  tr <- ss_trough(100, k, tau)
  withauc <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[1] * y[1], y[1]))
  }
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = tau,
                   ss = TRUE)
  got <- function(ext) {
    frm_ode(withauc, init = list(0, 0), times = 0, parms = list(k),
            events = ev, output = 1L, n_ss = 20L, ss_tol = Inf,
            ss_extrapolate = ext, atol = 1e-10, rtol = 1e-10)
  }
  e_trunc <- abs(got(FALSE) - tr) / tr
  e_ext <- abs(got(TRUE) - tr) / tr
  expect_lt(as.numeric(e_ext), as.numeric(e_trunc) / 1e4)
})

test_that("a state whose difference is exactly zero is left alone", {
  # An oral depot empties to the same trough every cycle once the
  # run-in has settled, so its successive differences are zero to the
  # last bit and an undamped ratio would be 0 / 0. Called directly
  # because that state is the one a caller never reads.
  ya <- c(0, 1)
  yb <- c(0, 2)
  yc <- c(0, 2.5)
  z <- ode_ss_extrapolate(ya, yb, yc, 1e-8, 1e-8)
  expect_true(all(is.finite(z$y)))
  expect_identical(z$y[1L], 0)
  # the state that IS moving still gets its tail: r = 0.5 here, so the
  # sum of the rest is 0.5 / (1 - 0.5) of the last difference
  expect_equal(z$y[2L], 2.5 + 0.5 * (0.5 / 0.5), tolerance = 1e-6)
})

test_that("a state that is growing between cycles is left alone", {
  # A ratio above 1 says the component is moving FURTHER each cycle, so
  # there is no tail to sum. Capping such a ratio at the top of its
  # range would read it as "nearly stopped", which is the worst reading
  # available: it multiplies the last difference by 1 / (1 - cap).
  z <- ode_ss_extrapolate(c(0, 1), c(0, 2), c(0, 4), 1e-8, 1e-8)
  expect_identical(z$y, c(0, 4))
  expect_gt(z$r[2L], 1)
})

test_that("the low gate takes tail away and never adds to it", {
  # With ya = 0, yb = 1/r and yc = 1/r + 1 the differences are 1/r and
  # 1, the damping is negligible and the ratio read is r, so what comes
  # back IS the correction factor. A gate may only remove tail: the
  # geometric tail is r / (1 - r) and the factor must not exceed it.
  # Shaped as S(x)/x rather than S(x) the gate overshot it by 1.2487x.
  fac <- function(r) {
    z <- ode_ss_extrapolate(0, 1 / r, 1 / r + 1, 1e-300, 1e-300)
    as.numeric(z$y) - (1 / r + 1)
  }
  rs <- seq(1e-4, ode_ss_rlow, length.out = 400)
  ratio <- vapply(rs, fac, 0) / (rs / (1 - rs))
  # a gate that saturates flat rises to the tail AT rlow and nowhere
  # exceeds it; the overshooting shape peaked inside the interval
  expect_identical(which.max(ratio), length(ratio))
  expect_lt(max(ratio), 1 + 1024 * .Machine$double.eps)
  expect_gt(ratio[length(ratio)], 0.99)
  # the junction is smooth: a corner holds the gap between the two
  # one-sided difference quotients up as the bracket shrinks
  gap <- function(eps) {
    f0 <- fac(ode_ss_rlow)
    abs((fac(ode_ss_rlow + eps) - f0) / eps -
          (f0 - fac(ode_ss_rlow - eps)) / eps)
  }
  expect_lt(gap(1e-5), gap(1e-3) / 10)
})

test_that("the guards cannot silence the warning they act on", {
  skip_if_not_installed("RTMBode")
  # Every guard works by driving the correction to exactly zero, so a
  # report built out of the correction alone is silent on exactly the
  # states the guards exist for. Three constructions, both arms: a state
  # with no steady state read AS THE OUTPUT, a run-in that oscillates,
  # and a run-in too short for the tail to be read at all.
  auc <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[2] * y[1], p[2] * y[1] - p[1] * y[2], y[2]))
  }
  osc <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(y[2], -p[1] * p[1] * y[1] - 2 * p[2] * p[1] * y[2]))
  }
  ev <- function(ii) data.frame(time = 0, state = 1L, value = 100,
                                ii = ii, ss = TRUE)
  for (ext in c(TRUE, FALSE)) {
    expect_warning(
      frm_ode(auc, init = list(0, 0, 0), times = c(0, 6),
              parms = list(0.05, 1), events = ev(12), output = 3L,
              n_ss = 20L, ss_extrapolate = ext, atol = 1e-10,
              rtol = 1e-10),
      "run-in", label = paste("AUC output, ss_extrapolate =", ext))
    expect_warning(
      frm_ode(osc, init = list(0, 0), times = c(0, 0.5),
              parms = list(3, 0.05), events = ev(1), output = 1L,
              n_ss = 20L, ss_extrapolate = ext, atol = 1e-10,
              rtol = 1e-10),
      "run-in", label = paste("oscillator, ss_extrapolate =", ext))
    expect_warning(
      frm_ode(decay1, init = list(0), times = 0, parms = list(0.01),
              events = ev(12), n_ss = 1L, ss_extrapolate = ext,
              atol = 1e-10, rtol = 1e-10),
      "run-in", label = paste("n_ss = 1, ss_extrapolate =", ext))
  }
})

test_that("a long run-in in the stand-down band still warns", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # `?frm_ode` sends a long-half-life user to a large `n_ss`, and that
  # is where the warning used to go quiet: the ratio sits ABOVE the
  # stand-down point, so part of the tail is deliberately not summed,
  # and successive extrapolants then agree with each other while both
  # sit short of the limit. The report adds the declined part back.
  k <- log(2) / 1653          # r = 0.98999, above ode_ss_rcap
  tau <- 24
  tr <- ss_trough(100, k, tau)
  ev <- data.frame(time = 0, value = 100, ii = tau, ss = TRUE)
  got <- function(n) {
    msg <- NULL
    v <- withCallingHandlers(
      frm_ode(decay1, init = list(0), times = 0, parms = list(k),
              events = ev, n_ss = n, atol = 1e-12, rtol = 1e-12),
      warning = function(w) {
        msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    list(err = abs(as.numeric(v) - tr) / tr, msg = msg)
  }
  z <- got(650L)
  # the error is real and above the default ss_tol, so silence here
  # would be a miss
  expect_gt(z$err, 1e-6)
  expect_true(!is.null(z$msg))
  said <- as.numeric(sub(".*about ([0-9.e+-]+) [(]relative.*", "\\1",
                         z$msg))
  # and what it says is the size of that error, not a token
  expect_gt(said / z$err, 0.5)
  expect_lt(said / z$err, 2)
  # far enough out it is right and says nothing; one call, because
  # every cycle is a solve and this test runs in every variant build
  far <- got(2000L)
  expect_lt(far$err, 1e-6)
  expect_null(far$msg)
})

test_that("ss_extrapolate = FALSE does not build the correction", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # Built outside the branch, the truncated arm puts the whole
  # correction on the tape and throws the value away. The node count is
  # load-independent, which a clock is not; the printed operation stack
  # is one line per node.
  ev <- data.frame(time = 0, value = 100, ii = 12, ss = TRUE)
  nodes <- function(ext) {
    tp <- RTMB::MakeTape(function(th) {
      "c" <- RTMB::ADoverload("c")
      sum(frm_ode(decay1, init = list(0), times = c(0, 6),
                  parms = list(exp(th[1])), events = ev, n_ss = 20L,
                  ss_tol = Inf, ss_extrapolate = ext))
    }, log(0.01))
    length(capture.output(tp$print(depth = 1))) - 1L
  }
  a <- nodes(FALSE)
  b <- nodes(TRUE)
  expect_gt(b, a)
  # and the difference is one correction, not a rounding difference
  expect_gt(b - a, 0.1 * a)
})

test_that("the tail correction leaves the objective differentiable", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # The correction stands down as the measured ratio approaches 1, and
  # a HARD cap there makes the objective C0 but not C1: the two
  # one-sided derivatives were -0.09 and +150 and did not converge as
  # the bracket tightened. The shape that replaced it is C1, so the gap
  # falls with the bracket.
  ke <- 0.15; k12 <- 0.3; ka <- 1.0; ii <- 24
  # k21 that puts exp(-lambda_z * ii) exactly at the stand-down point
  bigL <- -log(ode_ss_rcap) / ii
  k0 <- (bigL * bigL - (ke + k12) * bigL) / (bigL - ke)
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = ii,
                   ss = TRUE)
  dyn <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[4] * y[1], p[4] * y[1] - (p[1] + p[2]) * y[2] +
             p[3] * y[3], p[2] * y[2] - p[3] * y[3]))
  }
  g <- function(lk) {
    tp <- RTMB::MakeTape(function(th) {
      "c" <- RTMB::ADoverload("c")
      sum(frm_ode(dyn, init = list(0, 0, 0), times = c(0, 6, 23.9),
                  parms = list(ke, k12, exp(th[1]), ka), events = ev,
                  output = 2L, n_ss = 20L, ss_tol = Inf,
                  atol = 1e-12, rtol = 1e-12))
    }, lk)
    as.numeric(tp$jacfun()(lk))
  }
  gap <- function(eps) abs(g(log(k0 * (1 + eps))) -
                             g(log(k0 * (1 - eps))))
  wide <- gap(1e-3)
  tight <- gap(1e-4)
  # a derivative discontinuity holds the gap up as the bracket shrinks;
  # the bar is a ratio to what the same run measured
  expect_lt(tight, wide / 5)
})

test_that("the run-in warning reports a distance to the limit", {
  skip_if_not_installed("RTMBode")
  # The number the warning used to print was the movement between the
  # last two cycles, which is smaller than the distance to the limit by
  # about 1 / (k * ii) and so understated the error most where it was
  # largest. Here 1 / (k * ii) is about 8.
  k <- 0.01
  tau <- 12
  tr <- ss_trough(100, k, tau)
  ev <- data.frame(time = 0, value = 100, ii = tau, ss = TRUE)
  msg <- NULL
  got <- withCallingHandlers(
    frm_ode(decay1, init = list(0), times = 0, parms = list(k),
            events = ev, n_ss = 20L, ss_extrapolate = FALSE,
            atol = 1e-10, rtol = 1e-10),
    warning = function(w) {
      msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })
  expect_true(!is.null(msg))
  said <- as.numeric(sub(".*about ([0-9.e+-]+) [(]relative.*", "\\1",
                         msg))
  true <- abs(as.numeric(got) - tr) / tr
  # what it prints is the distance, so the ratio to the measured
  # distance is of order 1 rather than of order k * ii
  expect_gt(said / true, 0.5)
  expect_lt(said / true, 2)
  # and it is silent when the tail is summed, because then the answer
  # really is right
  expect_silent(
    frm_ode(decay1, init = list(0), times = 0, parms = list(k),
            events = ev, n_ss = 20L, atol = 1e-10, rtol = 1e-10)
  )
})

test_that("ss_extrapolate is validated", {
  skip_if_not_installed("RTMBode")
  expect_error(
    frm_ode(decay1, init = list(0), times = 0, parms = list(0.2),
            events = data.frame(time = 0, value = 100, ii = 12,
                                ss = TRUE), ss_extrapolate = "yes"),
    "`ss_extrapolate` must be TRUE or FALSE")
  expect_error(
    frm_ode(decay1, init = list(0), times = 0, parms = list(0.2),
            events = data.frame(time = 0, value = 100, ii = 12,
                                ss = TRUE), ss_extrapolate = NA),
    "`ss_extrapolate` must be TRUE or FALSE")
})

test_that("the gradient through a steady-state run-in is exact", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  obs <- c(0, 2, 6, 11)
  x0 <- c(log(0.2), 0.4)
  tp <- RTMB::MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    sum(frm_ode(decay1, init = list(0), times = obs,
                parms = list(exp(th[1])),
                events = data.frame(time = 0, value = 100, ii = 12,
                                    ss = TRUE),
                event_scale = 1 / (1 + exp(-th[2])),
                n_ss = 25, atol = 1e-12, rtol = 1e-12))
  }, x0)
  ref <- function(th) {
    k <- exp(th[1])
    D <- 100 / (1 + exp(-th[2]))
    tr <- ss_trough(D, k, 12)
    v <- (tr + D) * exp(-k * obs)
    v[1] <- tr
    sum(v)
  }
  expect_equal(tp(x0), ref(x0), tolerance = 1e-6)
  expect_equal(as.numeric(tp$jacfun()(x0)), central_fd(ref, x0),
               tolerance = 1e-5)
})

test_that("the tail correction fixes the gradient, not only the value", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  # The truncated run-in differentiates its own answer correctly and
  # that answer is the wrong one, which is what makes the shortfall a
  # wrong answer during a FIT rather than a slow one. k * ii = 0.12.
  tau <- 12
  obs <- c(0, 3, 9)
  x0 <- log(0.01)
  ev <- data.frame(time = 0, value = 100, ii = tau, ss = TRUE)
  ref <- function(th) {
    k <- exp(th[1])
    tr <- ss_trough(100, k, tau)
    v <- (tr + 100) * exp(-k * obs)
    v[1] <- tr
    sum(v)
  }
  g_ref <- central_fd(ref, x0)
  g <- function(ext) {
    tp <- RTMB::MakeTape(function(th) {
      "c" <- RTMB::ADoverload("c")
      sum(frm_ode(decay1, init = list(0), times = obs,
                  parms = list(exp(th[1])), events = ev, n_ss = 20L,
                  ss_tol = Inf, ss_extrapolate = ext, atol = 1e-10,
                  rtol = 1e-10))
    }, x0)
    as.numeric(tp$jacfun()(x0))
  }
  e_trunc <- abs(g(FALSE) - g_ref) / abs(g_ref)
  e_ext <- abs(g(TRUE) - g_ref) / abs(g_ref)
  expect_lt(e_ext, e_trunc / 1e3)
})

test_that("steady-state rows are validated", {
  skip_if_not_installed("RTMBode")
  call_ev <- function(ev) {
    frm_ode(decay1, init = list(0), times = c(1, 2), parms = list(0.2),
            events = ev)
  }
  expect_error(call_ev(data.frame(time = 0, value = 1, ss = TRUE)),
               "whose `ii` is not positive")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = 6, ss = TRUE,
                                  method = "replace")),
               "whose method is not")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = 6, ss = TRUE,
                                  duration = 9)),
               "longer than")
  expect_error(call_ev(data.frame(time = c(0, 3), value = 1, ii = 6,
                                  ss = TRUE)),
               "more than one row `ss = TRUE`")
  expect_error(call_ev(data.frame(time = 0, value = 1, ii = 6,
                                  ss = "yes")),
               "must be TRUE/FALSE")
})

test_that("a steady-state population fit recovers the truth", {
  skip_if_not_installed("RTMBode")
  skip_on_cran()
  set.seed(31)
  n_id <- 8
  tau <- 12
  tt <- c(0, 0.5, 2, 4, 8, 11.9)
  d <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                  time = rep(tt, n_id))
  k_i <- exp(stats::rnorm(n_id, log(0.2), 0.25))[as.integer(d$id)]
  tr <- ss_trough(100, k_i, tau)
  mu <- (tr + 100) * exp(-k_i * d$time)
  mu[d$time == 0] <- tr[d$time == 0]
  d$conc <- mu + stats::rnorm(nrow(d), 0, 0.5)

  fit <- frm(
    bf(conc ~ frm_ode(decay1, init = list(0), times = time,
                      parms = list(exp(lk)), group = id, output = 1L,
                      events = data.frame(time = 0, value = 100, ii = 12,
                                          ss = TRUE)),
       lk ~ 1 + (1 | id), nl = TRUE) + gaussian(),
    data = d, start = list(beta = log(0.25)))
  fx <- unlist(fixef(fit))
  expect_equal(unname(fx[["lk.(Intercept)"]]), log(0.2), tolerance = 0.15)
  expect_lt(unname(exp(fx[["sigma.(Intercept)"]])), 1.5)
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
