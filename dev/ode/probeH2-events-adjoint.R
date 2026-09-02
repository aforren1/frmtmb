# Probe H2: is the ADJOINT correct through deSolve events?
#
# This is the wrong-likelihood question. RTMBode::ode() forwards `...` to
# deSolve::ode() in both branches: once for the plain numeric solve, and
# again (inside ODEadjoint$updateSolution, and a second time for the
# augmented sensitivity system built by augment()) on the AD path.
#
# Structural prediction before running anything. The augmented state
# vector is c(y, dy...), so an event naming state k by INDEX still lands
# on the real state, but nothing jumps the sensitivity block:
#   method = "add"      y := y + a       d/dp (y + a) = dy/dp   -> correct
#   method = "replace"  y := v           d/dp v       = 0       -> WRONG
#   method = "multiply" y := f y         d/dp (f y)   = f dy/dp -> WRONG
# Sections 3 and 4 test that. Section 0 gets there first with a blunter
# answer.

suppressPackageStartupMessages({
  library(RTMB)
  library(RTMBode)
})
source("dev/ode/pk-common.R")

say <- function(...) cat(..., "\n", sep = "")
hr <- function(t) say("\n===== ", t, " =====")

D <- 100
dose_times <- c(0, 6, 12, 18)
ev_times <- dose_times[-1]
obs <- c(1, 3, 5, 7, 9, 11, 13, 17, 20, 24)
grid <- sort(unique(c(0, obs, ev_times)))
theta0 <- c(log(1.0), log(0.2), log(10))

ev_add <- data.frame(var = 1, time = ev_times, value = D, method = "add")

multi_dose_analytic <- function(t, ka, ke, V, D, dose_times) {
  vapply(t, function(ti) {
    u <- ti - dose_times[dose_times <= ti]
    if (!length(u)) return(0)
    sum(D * ka / (V * (ka - ke)) * (exp(-ke * u) - exp(-ka * u)))
  }, 0)
}

cfd <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}

report <- function(label, g_ad, g_ref) {
  say(label)
  say("  AD  : ", paste(format(g_ad, digits = 8), collapse = "  "))
  say("  ref : ", paste(format(g_ref, digits = 8), collapse = "  "))
  rel <- abs(g_ad - g_ref) / pmax(abs(g_ref), 1e-8)
  say("  max abs diff = ", format(max(abs(g_ad - g_ref)), digits = 4),
      "   max rel = ", format(max(rel), digits = 4))
  say("  VERDICT: ", if (max(rel) < 1e-4) "MATCHES" else "*** MISMATCH ***")
  invisible(max(rel))
}

f_ana <- function(th) {
  sum(multi_dose_analytic(grid[-1], exp(th[1]), exp(th[2]), exp(th[3]),
                          D, dose_times))
}

# --------------------------------------------------------------------
hr("0. does the AD path accept events = list(data = ) at all?")

solve_ad <- function(theta, events = NULL) {
  "c" <- ADoverload("c")
  p <- c(exp(theta[1]), exp(theta[2]), exp(theta[3]))
  args <- list(y = c(D, 0), times = grid, func = pk_dyn, parms = p,
               method = "lsoda", atol = 1e-10, rtol = 1e-10)
  if (!is.null(events)) args$events <- list(data = events)
  sol <- do.call(RTMBode::ode, args)
  sol[-1, 3]
}

tp <- tryCatch(MakeTape(function(th) sum(solve_ad(th, ev_add)), theta0),
               error = function(e) e)
if (inherits(tp, "error")) {
  say("*** FAIL taping with events = list(data = ): ",
      conditionMessage(tp))
} else {
  say("OK taped. value = ", format(tp(theta0), digits = 10))
}

# also: named states, in case the failure is about name lookup
tp_nm <- tryCatch(
  MakeTape(function(th) {
    "c" <- ADoverload("c")
    p <- c(exp(th[1]), exp(th[2]), exp(th[3]))
    sol <- RTMBode::ode(y = c(depot = D, central = 0), times = grid,
                        func = pk_dyn, parms = p, method = "lsoda",
                        atol = 1e-10, rtol = 1e-10,
                        events = list(data = data.frame(
                          var = "depot", time = ev_times, value = D,
                          method = "add")))
    sum(sol[-1, 3])
  }, theta0), error = function(e) e)
say("named-state spelling: ",
    if (inherits(tp_nm, "error")) paste("FAIL:", conditionMessage(tp_nm))
    else "OK")

# --------------------------------------------------------------------
hr("0b. why: the state vector RTMBode hands deSolve has no names")

ns <- asNamespace("RTMBode")
Ft <- ns$func2tape(pk_dyn, c(depot = D, central = 0), c(1, .2, 10))
say("names(F$par()) = ", if (is.null(names(Ft$par()))) "NULL"
    else paste(names(Ft$par()), collapse = ","))
Fi <- ns$addInfo(Ft, times = grid)
say("names(info$augstate) [order 0, the forward solve] = ",
    if (is.null(names(attr(Fi, "info")$augstate))) "NULL"
    else paste(names(attr(Fi, "info")$augstate), collapse = ","))
Gi <- ns$augment(Fi)
say("names(info$augstate) [order 1, the sensitivity solve] = ",
    paste(names(attr(Gi, "info")$augstate), collapse = ","))
say("deSolve::checkevents() bounds the event `var` index by ",
    "length(attr(y, 'names')). With NULL names that bound is 0, hence ",
    "'should be < 0'. This is an upstream RTMBode defect: ",
    "ODEadjoint()'s updateSolution() passes an unnamed augstate.")

# --------------------------------------------------------------------
hr("0c. workaround A: events = list(func = , time = )")

# checkevents() returns before it ever looks at the state names when an
# event FUNCTION is given, so this route sidesteps the missing names.
# The function is handed the AUGMENTED state, whose first nstate entries
# are the real states and whose tail is the sensitivity block. Adding a
# constant to a state leaves d state / d parm untouched, so touching only
# the first nstate entries is exactly right at every order.
ev_func <- function(t, y, p) {
  j <- which(abs(t - ev_times) < 1e-8)
  if (length(j)) y[1] <- y[1] + D
  y
}
tp_f <- tryCatch(
  MakeTape(function(th) {
    "c" <- ADoverload("c")
    p <- c(exp(th[1]), exp(th[2]), exp(th[3]))
    sol <- RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn,
                        parms = p, method = "lsoda",
                        atol = 1e-10, rtol = 1e-10,
                        events = list(func = ev_func, time = ev_times))
    sum(sol[-1, 3])
  }, theta0), error = function(e) e)
if (inherits(tp_f, "error")) {
  say("FAIL: ", conditionMessage(tp_f))
} else {
  say("taped. value = ", format(tp_f(theta0), digits = 10),
      "   analytic = ", format(f_ana(theta0), digits = 10))
  report("d/dtheta via events$func (ADD by hand)",
         as.numeric(tp_f$jacfun()(theta0)), cfd(f_ana, theta0))
}

# --------------------------------------------------------------------
hr("0d. workaround B: segmented solve, dose applied between segments")

# Split [t0, tmax] at the dose times and chain one ode() call per
# interval, carrying the end state forward and adding the dose in plain
# RTMB arithmetic. RTMBode is already differentiable through y0 (that is
# the whole point of the augmented system), so this is adjoint-correct by
# construction, with no event machinery involved at all.
seg_solve <- function(theta, doses_t, doses_a, y0, want) {
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  p <- c(exp(theta[1]), exp(theta[2]), exp(theta[3]))
  brk <- sort(unique(c(0, doses_t, max(want))))
  out <- numeric(length(want))
  y <- y0
  for (s in seq_len(length(brk) - 1L)) {
    a <- brk[s]; b <- brk[s + 1L]
    k <- which(abs(doses_t - a) < 1e-12)
    if (length(k)) y[1] <- y[1] + sum(doses_a[k])
    ts <- want[want > a & want <= b]
    gg <- c(a, ts, if (!length(ts) || ts[length(ts)] < b) b)
    sol <- RTMBode::ode(y = y, times = gg, func = pk_dyn, parms = p,
                        method = "lsoda", atol = 1e-10, rtol = 1e-10)
    if (length(ts)) out[match(ts, want)] <- sol[seq_along(ts) + 1L, 3]
    y <- c(sol[nrow(sol), 2], sol[nrow(sol), 3])
  }
  out
}
tp_s <- tryCatch(
  MakeTape(function(th)
    sum(seg_solve(th, dose_times, rep(D, 4), c(0, 0), obs)), theta0),
  error = function(e) e)
if (inherits(tp_s, "error")) {
  say("FAIL: ", conditionMessage(tp_s))
} else {
  v <- tp_s(theta0)
  a <- sum(multi_dose_analytic(obs, exp(theta0[1]), exp(theta0[2]),
                               exp(theta0[3]), D, dose_times))
  say("value = ", format(v, digits = 10), "  analytic = ",
      format(a, digits = 10), "  diff = ", format(abs(v - a), digits = 3))
  f_obs <- function(th) sum(multi_dose_analytic(
    obs, exp(th[1]), exp(th[2]), exp(th[3]), D, dose_times))
  report("d/dtheta via segmented solve",
         as.numeric(tp_s$jacfun()(theta0)), cfd(f_obs, theta0))
}

# --------------------------------------------------------------------
hr("0e. segmented solve: PARAMETER-DEPENDENT dose (bioavailability)")

if (!inherits(tp_s, "error")) {
  f_bio <- function(th) {
    Fb <- plogis(th[4])
    sum(seg_solve(th[1:3], dose_times, rep(Fb * D, 4), c(0, 0), obs))
  }
  tp_b <- tryCatch(MakeTape(f_bio, c(theta0, 0.5)), error = function(e) e)
  if (inherits(tp_b, "error")) {
    say("FAIL: ", conditionMessage(tp_b))
  } else {
    f_bio_num <- function(th)
      sum(multi_dose_analytic(obs, exp(th[1]), exp(th[2]), exp(th[3]),
                              plogis(th[4]) * D, dose_times))
    say("value = ", format(tp_b(c(theta0, 0.5)), digits = 10),
        "  analytic = ", format(f_bio_num(c(theta0, 0.5)), digits = 10))
    report("d/d(lka,lke,lV,logitF), dose = F * amt",
           as.numeric(tp_b$jacfun()(c(theta0, 0.5))),
           cfd(f_bio_num, c(theta0, 0.5)))
  }
}

# --------------------------------------------------------------------
hr("1. ADD events via events$data - blocked upstream, retest if fixed")
say(if (inherits(tp, "error")) "skipped (section 0 failed)" else "see 0")

# --------------------------------------------------------------------
hr("2. the sensitivity-block question, answered directly")

# Even with the names bug fixed, replace/multiply would still be
# suspect. Emulate what a *fixed* events$data would do by using
# events$func to touch ONLY the first nstate entries (which is what an
# indexed events$data row does) and compare against the truth.
mk <- function(meth, val) {
  function(t, y, p) {
    if (any(abs(t - ev_times) < 1e-8)) {
      y[1] <- switch(meth, replace = val, multiply = val * y[1],
                     add = y[1] + val)
    }
    y
  }
}
num_ref <- function(th, meth, val) {
  p <- c(exp(th[1]), exp(th[2]), exp(th[3]))
  sol <- RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn, parms = p,
                      method = "lsoda", atol = 1e-10, rtol = 1e-10,
                      events = list(data = data.frame(
                        var = 1, time = ev_times, value = val,
                        method = meth)))
  sum(sol[-1, 3])
}
for (m in c("add", "replace", "multiply")) {
  val <- switch(m, add = D, replace = 50, multiply = 2)
  tt <- tryCatch(
    MakeTape(function(th) {
      "c" <- ADoverload("c")
      p <- c(exp(th[1]), exp(th[2]), exp(th[3]))
      sol <- RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn,
                          parms = p, method = "lsoda",
                          atol = 1e-10, rtol = 1e-10,
                          events = list(func = mk(m, val),
                                        time = ev_times))
      sum(sol[-1, 3])
    }, theta0), error = function(e) e)
  if (inherits(tt, "error")) {
    say(m, ": taping FAILED: ", conditionMessage(tt))
    next
  }
  say("\n-- method = ", m, " --")
  say("  value AD = ", format(tt(theta0), digits = 10),
      "   value numeric events$data = ",
      format(num_ref(theta0, m, val), digits = 10))
  report(paste0("  d/dtheta, state-only jump, ", m),
         as.numeric(tt$jacfun()(theta0)),
         cfd(function(th) num_ref(th, m, val), theta0))
}

say("\ndone.")
