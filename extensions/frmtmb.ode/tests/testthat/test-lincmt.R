# frm_lincmt(): the analytic linear compartment model.
#
# The reference is frm_ode() on the same schedule at the same
# parameters, which is an independent implementation of the same model
# through a numerical solve. No tolerance below is a constant. The
# solver is run TWICE, at `tol` and at `tol / 1000`, and the difference
# between those two runs is what the solver's own tolerance costs on
# that schedule; the closed form is then required to agree with the
# tighter run by a small multiple of that, floored at a few ulp of the
# trajectory's own scale. So a schedule the solver finds hard gets a
# looser bar and a schedule it finds easy gets a tighter one, both
# measured by the run.
#
# The one place a constant would be unavoidable is the analytic limit
# `ka == ke`, and there the reference is exact arithmetic rather than a
# solve, so the bar is ulp of the value.

# every reference system below is written in AMOUNTS, so a
# concentration is the frm_ode() result divided by V
pk1_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
}
pk2_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[4] * y[1],
         p[4] * y[1] - (p[1] + p[2]) * y[2] + p[3] * y[3],
         p[2] * y[2] - p[3] * y[3]))
}
iv1_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}

# The bar: what the solver's own tolerance costs on THIS schedule,
# with a floor at a few ulp of the trajectory's scale.
solver_bar <- function(loose, tight, k = 5) {
  spread <- max(abs(loose - tight))
  k * max(spread, 8 * .Machine$double.eps * max(abs(tight)))
}

fd_grad <- function(f, x, h) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[[j]] <- xp[[j]] + h
    xm <- x; xm[[j]] <- xm[[j]] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}

tt <- c(0, 0.5, 1, 2, 4, 8, 12, 18, 24, 30, 36, 48)

# --- identity with frm_ode() across the schedule space ---------------

test_that("one compartment with a depot matches frm_ode() on every
           schedule the closed form covers", {
  skip_if_not_installed("RTMBode")
  P <- list(ka = 1.1, ke = 0.2, V = 10)
  sched <- list(
    "init bolus" = list(init = list(depot = 100), events = NULL),
    "single dose" = list(
      events = data.frame(time = 2, state = "depot", value = 100)),
    "addl and ii" = list(
      events = data.frame(time = 0, state = "depot", value = 100,
                          ii = 8, addl = 5L)),
    "steady state" = list(
      events = data.frame(time = 0, state = "depot", value = 100,
                          ii = 8, ss = TRUE)),
    "steady state plus addl" = list(
      events = data.frame(time = c(0, 8), state = "depot", value = 100,
                          ii = c(8, 8), addl = c(0L, 5L),
                          ss = c(TRUE, FALSE))),
    "infusion into the central compartment" = list(
      events = data.frame(time = 1, state = "central", value = 100,
                          duration = 3)),
    "repeated infusion" = list(
      events = data.frame(time = 0, state = "central", value = 100,
                          duration = 3, ii = 8, addl = 5L)),
    "infusion at steady state" = list(
      events = data.frame(time = 0, state = "central", value = 100,
                          duration = 3, ii = 8, ss = TRUE)),
    "reset to zero" = list(
      init = list(depot = 100),
      events = data.frame(time = c(12, 12), state = c(NA, "depot"),
                          value = c(0, 100),
                          method = c("reset", "add"))))
  for (nm in names(sched)) {
    s <- sched[[nm]]
    ini <- s[["init"]]
    i0 <- list(if (is.null(ini)) 0 else ini[["depot"]], 0)
    a <- frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
                    init = ini, events = s[["events"]], n_ss = 20L)
    ode_at <- function(tol) {
      frm_ode(pk1_dyn, init = i0, times = tt,
              parms = list(P$ka, P$ke),
              states = c("depot", "central"), output = "central",
              events = s[["events"]], n_ss = 20L,
              atol = tol, rtol = tol) / P$V
    }
    loose <- ode_at(1e-9)
    tight <- ode_at(1e-12)
    expect_lt(max(abs(a - tight)), solver_bar(loose, tight),
              label = paste("1 cmt, depot,", nm))
  }
})

test_that("two and three compartments match frm_ode(), with and
           without a depot", {
  skip_if_not_installed("RTMBode")
  P <- list(ka = 1.1, ke = 0.2, k12 = 0.4, k21 = 0.1, V = 10)
  ev <- data.frame(time = c(0, 8), state = "depot", value = 100,
                   ii = c(8, 8), addl = c(0L, 3L), ss = c(TRUE, FALSE))
  # The DEFAULT arms are compared, and they agree because both reach
  # the limit: frm_lincmt() sums the geometric series and frm_ode()
  # sums the tail its run-in would otherwise truncate. Before frm_ode()
  # summed that tail this comparison was 3.5e-02 apart on this
  # schedule, and the test matched the two truncations to each other.
  a <- frm_lincmt(parms = P, times = tt, ncmt = 2, depot = TRUE,
                  events = ev)
  ode_at <- function(tol, ...) suppressWarnings({
    frm_ode(pk2_dyn, init = list(0, 0, 0), times = tt,
            parms = list(P$ke, P$k12, P$k21, P$ka),
            states = c("depot", "central", "peripheral1"),
            output = "central", events = ev, n_ss = 20L,
            atol = tol, rtol = tol, ...) / P$V
  })
  loose <- ode_at(1e-9); tight <- ode_at(1e-12)
  expect_lt(max(abs(a - tight)), solver_bar(loose, tight))
  # and the matched truncations still agree with each other, which is
  # what ss_extrapolate = FALSE is for
  b <- frm_lincmt(parms = P, times = tt, ncmt = 2, depot = TRUE,
                  events = ev, n_ss = 20L)
  loose2 <- ode_at(1e-9, ss_extrapolate = FALSE)
  tight2 <- ode_at(1e-12, ss_extrapolate = FALSE)
  expect_lt(max(abs(b - tight2)), solver_bar(loose2, tight2))
  # the two truncations are far from the limit, which is the defect
  expect_gt(max(abs(b - a)), 1e-3 * max(abs(a)))

  # three compartments, no depot, an intravenous infusion
  P3 <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05, k31 = 0.01,
             V = 10)
  d3 <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-(p[1] + p[2] + p[4]) * y[1] + p[3] * y[2] + p[5] * y[3],
           p[2] * y[1] - p[3] * y[2],
           p[4] * y[1] - p[5] * y[3]))
  }
  ev3 <- data.frame(time = 0, state = "central", value = 100,
                    duration = 2, ii = 8, addl = 4L)
  a3 <- frm_lincmt(parms = P3, times = tt, ncmt = 3, depot = FALSE,
                   events = ev3)
  ode3 <- function(tol) {
    frm_ode(d3, init = list(0, 0, 0), times = tt,
            parms = list(P3$ke, P3$k12, P3$k21, P3$k13, P3$k31),
            states = c("central", "peripheral1", "peripheral2"),
            output = "central", events = ev3, atol = tol,
            rtol = tol) / P3$V
  }
  l3 <- ode3(1e-9); t3 <- ode3(1e-12)
  expect_lt(max(abs(a3 - t3)), solver_bar(l3, t3))
})

test_that("groups, a per-group schedule and an estimated dose scale
           match frm_ode()", {
  skip_if_not_installed("RTMBode")
  d <- data.frame(id = factor(rep(1:3, each = 6)),
                  time = rep(c(0, 1, 4, 8, 12, 24), 3))
  sc <- c(0.6, 0.9, 1.2)[as.integer(d$id)]
  ev <- data.frame(group = as.character(rep(1:3, each = 2)),
                   time = rep(c(0, 8), 3), state = "depot",
                   value = rep(c(100, 50, 200), each = 2), ii = 8,
                   addl = rep(c(1L, 2L, 0L), each = 2))
  P <- list(ka = 1.1, ke = 0.2, k12 = 0.4, k21 = 0.1, V = 10)
  a <- frm_lincmt(parms = P, times = d$time, group = d$id, ncmt = 2,
                  depot = TRUE, events = ev, event_scale = sc)
  ode_at <- function(tol) {
    frm_ode(pk2_dyn, init = list(0, 0, 0), times = d$time,
            group = d$id, parms = list(P$ke, P$k12, P$k21, P$ka),
            states = c("depot", "central", "peripheral1"),
            output = "central", events = ev, event_scale = sc,
            atol = tol, rtol = tol) / P$V
  }
  loose <- ode_at(1e-9); tight <- ode_at(1e-12)
  expect_lt(max(abs(a - tight)), solver_bar(loose, tight))
})

test_that("the depot amount and the clearance spelling agree with what
           they stand for", {
  skip_if_not_installed("RTMBode")
  P <- list(ka = 1.1, ke = 0.2, V = 10)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                   addl = 3L)
  a <- frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
                  events = ev, output = "depot")
  ode_at <- function(tol) {
    frm_ode(pk1_dyn, init = list(0, 0), times = tt,
            parms = list(P$ka, P$ke),
            states = c("depot", "central"), output = "depot",
            events = ev, atol = tol, rtol = tol)
  }
  loose <- ode_at(1e-9); tight <- ode_at(1e-12)
  expect_lt(max(abs(a - tight)), solver_bar(loose, tight))

  # CL / V / Q2 / V2 is the same model written the other way
  cl <- list(CL = 0.2 * 10, V = 10, Q2 = 0.4 * 10, V2 = 0.4 * 10 / 0.1,
             ka = 1.1)
  rate <- list(ke = 0.2, V = 10, k12 = 0.4, k21 = 0.1, ka = 1.1)
  x <- frm_lincmt(parms = cl, times = tt, ncmt = 2, depot = TRUE,
                  events = ev)
  y <- frm_lincmt(parms = rate, times = tt, ncmt = 2, depot = TRUE,
                  events = ev)
  expect_lt(max(abs(x - y)), 16 * .Machine$double.eps * max(abs(y)))
})

# --- the numerics the item exists for --------------------------------

test_that("ka == ke is the analytic limit and not a NaN", {
  k <- 0.2
  V <- 10
  D <- 100
  a <- frm_lincmt(parms = list(ka = k, ke = k, V = V), times = tt,
                  ncmt = 1, depot = TRUE, init = list(depot = D))
  # the limit of D ka (exp(-ke t) - exp(-ka t)) / (V (ka - ke))
  want <- D * k * tt * exp(-k * tt) / V
  expect_true(all(is.finite(a)))
  expect_lt(max(abs(a - want)), 8 * .Machine$double.eps * max(want))
})

test_that("the closed form keeps its digits where the textbook form
           loses them, as ka approaches ke", {
  k <- 0.2; V <- 10; D <- 100
  u <- 6
  # the textbook difference of exponentials, for comparison
  naive <- function(ka) {
    D * ka / (V * (ka - k)) * (exp(-k * u) - exp(-ka * u))
  }
  want <- D * k * u * exp(-k * u) / V
  err <- function(f, d) abs(f(k + d) - want) / want
  ours <- function(ka) {
    frm_lincmt(parms = list(ka = ka, ke = k, V = V), times = u,
               ncmt = 1, depot = TRUE, init = list(depot = D))
  }
  # at a separation of 1e-12 the naive form has lost most of its digits
  # and the closed form has lost none; the bar is the naive form's own
  # measured error rather than a constant
  d <- 1e-12
  expect_gt(err(naive, d), 1e3 * err(ours, d))
  # and at zero separation the naive form is not a number at all
  expect_true(is.nan(naive(k)))
  expect_true(is.finite(ours(k)))
})

test_that("the gradient is finite and right where two rate constants
           meet, which is where the first version was not", {
  # This pins a defect that was in this file's own first version: the
  # value at ka == ke was correct but the tape returned 1.2e286 for
  # d/d(log ka), because (1 - exp(-z)) / z at z = 1e-300 differentiates
  # by the chain rule into two numbers of size 1e300 whose difference
  # is -1/2.
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                   addl = 3L)
  f <- function(th) {
    sum(frm_lincmt(parms = list(ka = exp(th[1]), ke = exp(th[2]),
                                V = exp(th[3])),
                   times = tt, ncmt = 1, depot = TRUE, events = ev))
  }
  th <- c(log(0.2), log(0.2), log(10))
  g <- as.numeric(RTMB::MakeTape(f, th)$jacobian(th))
  expect_true(all(is.finite(g)))
  # the gradient just off the degenerate point is what the gradient AT
  # it must look like: a run that reports 1e286 fails here by 280
  # orders of magnitude
  near <- c(log(0.2 * (1 + 1e-4)), log(0.2), log(10))
  gn <- as.numeric(RTMB::MakeTape(f, near)$jacobian(near))
  expect_lt(max(abs(g)), 2 * max(abs(gn)))
  # and it agrees with a central difference of the same function, to a
  # bar set by the difference's own step-size sensitivity
  g1 <- fd_grad(f, th, 1e-5)
  g2 <- fd_grad(f, th, 1e-4)
  expect_lt(max(abs(g - g1)),
            20 * max(max(abs(g1 - g2)),
                     .Machine$double.eps * max(abs(g1))))
})

test_that("a genuine eigenvalue collision keeps its value, and its
           gradient once the roots are resolved", {
  skip_if_not_installed("RTMBode")
  # Setting every rate constant equal does NOT collide the disposition
  # eigenvalues of a three-compartment model: at ke = k12 = k21 = k13 =
  # k31 = r the characteristic polynomial factors as
  # (L - r)(L^2 - 4 r L + r^2) and the roots are r, (2 - sqrt 3) r and
  # (2 + sqrt 3) r, which are distinct. What collides there is `ka`
  # with one eigenvalue, which is lincmt_diff()'s singularity and not
  # lincmt_disp()'s. The disposition collision needs constructing: k13
  # at 1e-300 decouples the third compartment, and k31 on the slow root
  # of the reduced quadratic makes that root double.
  #
  # lincmt_disp() reaches that through acos(), whose derivative is
  # infinite where a double root puts its argument at one, so the
  # gradient is NaN while the two roots are unresolved and finite once
  # they are. dev/lincmt/lincmt-f6.R measures the line: NaN while the
  # split has underflowed to 2.8e-16, finite and agreeing with a
  # central difference to 2.4e-09 from a split of 5.2e-09 upward. The
  # value is right on both sides. ?frm_lincmt's "Boundary" says so.
  ke <- 0.2; k12 <- 0.4; k21 <- 0.1; k13 <- 1e-300
  b <- ke + k12 + k21
  slow <- (b - sqrt(b * b - 4 * ke * k21)) / 2
  d <- frmtmb.ode:::lincmt_disp(3L, ke, k12, k21, k13, slow)
  lam <- sort(unlist(d[["lam"]]))
  # the collision is real: two of the three roots agree to rounding
  expect_lt(min(diff(lam)), 64 * .Machine$double.eps * max(lam))

  ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                   addl = 3L)
  d3 <- function(t, y, p) {
    "c" <- RTMB::ADoverload("c")
    list(c(-p[6] * y[1],
           p[6] * y[1] - (p[1] + p[2] + p[4]) * y[2] + p[3] * y[3] +
             p[5] * y[4],
           p[2] * y[2] - p[3] * y[3],
           p[4] * y[2] - p[5] * y[4]))
  }
  # the VALUE at the exact tangency, and just off it
  for (k31 in c(slow, slow * (1 + 1e-8))) {
    a <- frm_lincmt(parms = list(ke = ke, k12 = k12, k21 = k21,
                                 k13 = k13, k31 = k31, ka = 1.1,
                                 V = 10),
                    times = tt, ncmt = 3, depot = TRUE, events = ev)
    ode_at <- function(tol) {
      frm_ode(d3, init = list(0, 0, 0, 0), times = tt,
              parms = list(ke, k12, k21, k13, k31, 1.1),
              states = c("depot", "central", "peripheral1",
                         "peripheral2"),
              output = "central", events = ev, atol = tol,
              rtol = tol) / 10
    }
    loose <- ode_at(1e-9); tight <- ode_at(1e-12)
    expect_lt(max(abs(a - tight)), solver_bar(loose, tight))
  }

  # the GRADIENT once the two roots are resolved
  k31 <- slow * (1 + 1e-8)
  f <- function(th) {
    sum(frm_lincmt(parms = list(ke = exp(th[1]), k12 = k12, k21 = k21,
                                k13 = k13, k31 = exp(th[2]),
                                ka = exp(th[3]), V = 10),
                   times = tt, ncmt = 3, depot = TRUE, events = ev))
  }
  th <- c(log(ke), log(k31), log(1.1))
  gr <- as.numeric(RTMB::MakeTape(f, th)$jacobian(th))
  expect_true(all(is.finite(gr)))
  g1 <- fd_grad(f, th, 1e-6)
  g2 <- fd_grad(f, th, 1e-5)
  expect_lt(max(abs(gr - g1)),
            50 * max(max(abs(g1 - g2)),
                     .Machine$double.eps * max(abs(g1))))
})

test_that("the tape agrees with frm_ode()'s adjoint gradient", {
  skip_if_not_installed("RTMBode")
  ev <- data.frame(time = c(0, 12), state = "depot", value = 100,
                   ii = c(12, 12), addl = c(0L, 3L),
                   ss = c(TRUE, FALSE))
  fl <- function(th) {
    sum(frm_lincmt(parms = list(ka = exp(th[1]), ke = exp(th[2]),
                                V = exp(th[3])),
                   times = tt, ncmt = 1, depot = TRUE, events = ev,
                   n_ss = 20L))
  }
  fo <- function(tol) function(th) {
    sum(frm_ode(pk1_dyn, init = list(0, 0), times = tt,
                parms = list(exp(th[1]), exp(th[2])),
                states = c("depot", "central"), output = "central",
                events = ev, n_ss = 20L, atol = tol,
                rtol = tol)) / exp(th[3])
  }
  jac <- function(f, th) as.numeric(RTMB::MakeTape(f, th)$jacobian(th))
  for (th in list(c(log(1), log(0.2), log(10)),
                  c(log(0.2), log(0.2), log(10)))) {
    gl <- jac(fl, th)
    loose <- jac(fo(1e-9), th)
    tight <- jac(fo(1e-12), th)
    # the ADJOINT gradient carries the solve's own tolerance, so the
    # bar is what loosening that tolerance by a thousand costs it
    expect_lt(max(abs(gl - tight)), solver_bar(loose, tight))
  }
})

test_that("a finite n_ss converges on the exact limit geometrically", {
  P <- list(ka = 1.0, ke = 0.15, V = 20)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                   ss = TRUE)
  ts <- 0.5
  exact <- frm_lincmt(parms = P, times = ts, ncmt = 1, depot = TRUE,
                      events = ev)
  gap <- vapply(c(5L, 10L, 15L), function(n) {
    abs(frm_lincmt(parms = P, times = ts, ncmt = 1, depot = TRUE,
                   events = ev, n_ss = n) - exact)
  }, 0)
  # the shortfall after n cycles is the accumulation factor to the
  # power n, so each five more cycles multiplies it by exp(-5 ke ii)
  ratio <- exp(-5 * P$ke * 12)
  expect_lt(gap[[2L]], 2 * ratio * gap[[1L]])
  expect_lt(gap[[3L]], 2 * ratio * gap[[2L]])
  # and at the default it is below the trajectory's own rounding
  expect_lt(abs(frm_lincmt(parms = P, times = ts, ncmt = 1,
                           depot = TRUE, events = ev, n_ss = 30L) -
                  exact),
            8 * .Machine$double.eps * exact)
})

test_that("a drug that is not eliminated accumulates without bound and
           without a NaN", {
  # `n_ss = Inf` is a limit, and the limit of an unbounded series is
  # unbounded. What must not happen is Inf or NaN, because an optimizer
  # backs out of a large finite objective and cannot back out of a NaN.
  ev <- data.frame(time = 0, state = "central", value = 100, ii = 12,
                   ss = TRUE)
  small <- vapply(c(1e-6, 1e-12, 1e-200, 0), function(k) {
    frm_lincmt(parms = list(ke = k, V = 1), times = 3, ncmt = 1,
               depot = FALSE, events = ev)
  }, 0)
  expect_true(all(is.finite(small)))
  # each factor of a million off the elimination rate is a factor of a
  # million more accumulation, until the offset in the primitive takes
  # over: the series is 1 / (1 - exp(-ke ii)), so the ratio is the
  # ratio of the rates
  expect_gt(small[[2L]] / small[[1L]], 1e5)
  # and the gradient is finite there too
  g <- RTMB::MakeTape(function(th) {
    frm_lincmt(parms = list(ke = exp(th[1]), V = 1), times = 3,
               ncmt = 1, depot = FALSE, events = ev)
  }, log(1e-6))$jacobian(log(1e-6))
  expect_true(all(is.finite(as.numeric(g))))
})

test_that("an observation at a dose time reads the trough, and one at
           t0 reads init", {
  P <- list(ka = 1.1, ke = 0.2, V = 10)
  ev <- data.frame(time = c(0, 8), state = "depot", value = 100)
  a <- frm_lincmt(parms = P, times = c(0, 8), ncmt = 1, depot = TRUE,
                  init = list(central = 30), events = ev,
                  output = "central")
  # at t0 the dose at t0 has not been given: the reading is init
  expect_identical(a[[1L]], 30)
  # at t = 8 the second dose has not been given either, so the reading
  # is what the first dose alone produced
  b <- frm_lincmt(parms = P, times = 8, ncmt = 1, depot = TRUE,
                  init = list(central = 30),
                  events = data.frame(time = 0, state = "depot",
                                      value = 100),
                  output = "central")
  expect_identical(a[[2L]], b[[1L]])
})

# --- what it refuses -------------------------------------------------

test_that("the parts of the grammar superposition cannot carry are
           refused by name", {
  P <- list(ka = 1.1, ke = 0.2, V = 10)
  call_ev <- function(ev, ...) {
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               events = ev, ...)
  }
  expect_error(
    call_ev(data.frame(time = 2, state = "depot", value = 100,
                       method = "replace")),
    "replace")
  expect_error(
    call_ev(data.frame(time = 2, state = "depot", value = 2,
                       method = "multiply")),
    "multiply")
  expect_error(
    call_ev(data.frame(time = 2, state = NA, value = 5,
                       method = "reset")),
    "reset")
  expect_error(
    call_ev(data.frame(time = 2, state = "depot", value = 100,
                       duration = 1)),
    "zero-order absorption")
  expect_error(
    frm_lincmt(parms = list(ka = 1, ke = 0.2, k12 = 0.4, k21 = 0.1,
                            V = 10),
               times = tt, ncmt = 2, depot = TRUE,
               events = data.frame(time = 0, state = "peripheral1",
                                   value = 100)),
    "peripheral compartment")
  expect_error(
    frm_lincmt(parms = list(ka = 1, ke = 0.2, k12 = 0.4, k21 = 0.1,
                            V = 10),
               times = tt, ncmt = 2, depot = TRUE,
               init = list(peripheral1 = 100)),
    "peripheral")
  expect_error(
    frm_lincmt(parms = P, times = tt, ncmt = 1, depot = TRUE,
               tv = list(0.2)),
    "`tv` is not available")
  # an infusion that spans a restart loses what it had still to deliver
  expect_error(
    call_ev(data.frame(time = c(0, 4), state = c("central", NA),
                       value = c(100, 0), duration = c(8, 0),
                       method = c("add", "reset"))),
    "infusion running at")
})

test_that("the parameter list is read strictly", {
  expect_error(frm_lincmt(parms = list(1, 0.2, 10), times = tt,
                          ncmt = 1, depot = TRUE),
               "must be named")
  expect_error(frm_lincmt(parms = list(ka = 1, ke = 0.2), times = tt,
                          ncmt = 1, depot = TRUE),
               "missing V")
  expect_error(frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10,
                                       k12 = 0.4), times = tt,
                          ncmt = 1, depot = TRUE),
               "k12")
  expect_error(frm_lincmt(parms = list(ka = 1, CL = 2, ke = 0.2,
                                       V = 10), times = tt, ncmt = 1,
                          depot = TRUE),
               "mixes the two parameterizations")
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10), times = tt,
                          ncmt = 4),
               "`ncmt` must be 1, 2 or 3")
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10), times = tt,
                          output = "depot"),
               "there is no depot")
  # frm_ode()'s per-row `output` reads a different state on each row
  # over one shared solve; this is the drop-in gap, named
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10), times = tt,
                          output = rep("central", length(tt))),
               "returns ONE column")
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10), times = tt,
                          n_ss = 0),
               "`n_ss` must be Inf")
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10), times = tt,
                          event_scale = 2),
               "nothing to scale")
  expect_error(frm_lincmt(parms = list(ke = 0.2, V = 10),
                          times = c(1, 2), t0 = 3),
               "before t0")
})

test_that("frm_lincmt() needs no solver backend", {
  # the claim that a closed-form fit has no RTMBode and no deSolve in
  # it, checked on the code rather than on the environment
  src <- unlist(lapply(
    c("frm_lincmt", "lincmt_resp", "lincmt_disp", "lincmt_phi",
      "lincmt_diff", "lincmt_e2", "lincmt_e2_ss", "lincmt_geo",
      "lincmt_rates", "lincmt_check_events", "lincmt_col",
      "lincmt_needed", "lincmt_relu"),
    function(nm) deparse(get(nm, envir = asNamespace("frmtmb.ode")))))
  expect_false(any(grepl("RTMBode|deSolve|ode_require_backend", src)))
})

# --- inside a formula ------------------------------------------------

test_that("a hierarchical fit through frm_lincmt() recovers its
           simulator", {
  skip_on_cran()
  set.seed(4321)
  ns <- 25L
  ts <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  lka <- log(1) + stats::rnorm(ns, 0, 0.3)
  lke <- log(0.2) + stats::rnorm(ns, 0, 0.25)
  d <- data.frame(id = factor(rep(seq_len(ns), each = length(ts))),
                  time = rep(ts, ns), dose = 100)
  i <- as.integer(d$id)
  ka <- exp(lka[i]); ke <- exp(lke[i])
  d$conc <- 100 * ka / (10 * (ka - ke)) *
    (exp(-ke * d$time) - exp(-ka * d$time)) +
    stats::rnorm(nrow(d), 0, 0.2)
  fit <- frm(
    bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                      V = exp(lV)),
                         times = time, group = id, ncmt = 1,
                         depot = TRUE, init = list(depot = dose)),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = d, start = list(beta = c(0, log(0.25), log(8))), se = TRUE)
  ci <- suppressWarnings(stats::confint(fit))
  want <- c("lka_(Intercept)" = log(1), "lke_(Intercept)" = log(0.2),
            "lV_(Intercept)" = log(10))
  for (nm in names(want)) {
    j <- which(rownames(ci) == nm)
    expect_length(j, 1L)
    expect_true(ci[j, 1L] <= want[[nm]] && ci[j, 2L] >= want[[nm]],
                label = nm)
  }
})

test_that("a dynamics input that varies inside a group is refused at
           frame assembly, for frm_lincmt() as for frm_ode()", {
  d <- data.frame(id = factor(rep(1:4, each = 4)),
                  time = rep(c(0.5, 1, 2, 4), 4), dose = 100,
                  conc = 1, phase = factor(rep(c("a", "a", "b", "b"),
                                               4)))
  expect_error(
    frm(bf(conc ~ frm_lincmt(parms = list(ka = exp(lka),
                                          ke = exp(lke),
                                          V = exp(lV)),
                             times = time, group = id, ncmt = 1,
                             depot = TRUE, init = list(depot = dose)),
           lka ~ 1, lke ~ 1 + phase, lV ~ 1, nl = TRUE) + gaussian(),
        data = d, start = list(beta = c(0, 0, 0, 0))),
    "not constant within 'id'")
})
