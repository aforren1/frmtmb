# Probe H3: time-dependent input. Infusions, forcings, and the trap.
#
# Task item (d) asked whether "branching on time is legal on the tape,
# since time is not a parameter". For plain deSolve it is. For RTMBode it
# is NOT, and the reason is visible in func2tape():
#
#   x <- numeric(1 + nstate + length(parms))    # ALL ZERO
#   MakeTape(function(typ) { t <- typ[1]; ... func(t, y, p) }, x)
#
# The derivative function is taped ONCE, at t = 0. `t` is a tape input,
# so arithmetic in t is differentiated correctly, but any R-level `if`,
# `which`, `approxfun` or other control flow on `t` is resolved at tape
# construction time, with t = 0, and frozen there for the whole
# integration. No error, no warning: a different model, silently.
#
# What actually happens is better than "frozen": `t` is an advector on
# the inner tape, RTMB refuses comparison on AD types, and the branch
# raises "Comparison is generally unsafe for AD types". Loudly, and on
# the NUMERIC path too, because func2tape() tapes the derivative function
# in both branches of RTMBode::ode(). Sections 1 and 2 pin that down.
# Section 4 gives the route that does work: a constant infusion rate
# carried as an extra PARAMETER, held constant over a segment of a
# segmented solve.
#
# Reference model, one-compartment IV infusion of rate R over [0, Tinf],
# elimination ke, volume V, amount A:
#   dA/dt = R(t) - ke A,  A(0) = 0
#   A(t) = R/ke (1 - exp(-ke t))                        t <= Tinf
#   A(t) = R/ke (1 - exp(-ke Tinf)) exp(-ke (t - Tinf))  t >  Tinf

suppressPackageStartupMessages({
  library(RTMB)
  library(RTMBode)
})

say <- function(...) cat(..., "\n", sep = "")
hr <- function(t) say("\n===== ", t, " =====")

R_inf <- 20; Tinf <- 4; ke0 <- 0.2
obs <- c(0.5, 1, 2, 3, 4, 5, 6, 8, 12)

inf_analytic <- function(t, R, ke, Tinf) {
  ifelse(t <= Tinf,
         R / ke * (1 - exp(-ke * t)),
         R / ke * (1 - exp(-ke * Tinf)) * exp(-ke * (t - Tinf)))
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
  say("  max rel = ", format(max(rel), digits = 4), "   VERDICT: ",
      if (max(rel) < 1e-4) "MATCHES" else "*** MISMATCH ***")
  invisible(max(rel))
}

# --------------------------------------------------------------------
hr("1. `if (t < Tinf)` inside dynamics: numeric vs taped")

dyn_if <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  rate <- if (t < p[2]) p[3] else 0
  list(c(rate - p[1] * y[1]))
}

ana <- inf_analytic(obs, R_inf, ke0, Tinf)
never_off <- R_inf / ke0 * (1 - exp(-ke0 * obs))   # rate on forever

# plain deSolve, no RTMBode: the branch is live and the answer is right
num_de <- deSolve::ode(y = c(0), times = c(0, obs), func = dyn_if,
                       parms = c(ke0, Tinf, R_inf), method = "lsoda",
                       atol = 1e-10, rtol = 1e-10)[-1, 2]
say("plain deSolve  max |ode - analytic| = ",
    format(max(abs(num_de - ana)), digits = 4), "  (works)")

# RTMBode NUMERIC path: func2tape() still tapes the derivative function,
# so `t` is an advector even here
num <- tryCatch(
  RTMBode::ode(y = c(0), times = c(0, obs), func = dyn_if,
               parms = c(ke0, Tinf, R_inf), method = "lsoda",
               atol = 1e-10, rtol = 1e-10)[-1, 2],
  error = function(e) e)
say("RTMBode numeric: ", if (inherits(num, "error"))
  paste("ERROR -", conditionMessage(num)) else "ok")

# RTMBode AD path
ad <- tryCatch(
  MakeTape(function(th) {
    "c" <- ADoverload("c")
    sol <- RTMBode::ode(y = c(0), times = c(0, obs), func = dyn_if,
                        parms = c(exp(th[1]), Tinf, R_inf),
                        method = "lsoda", atol = 1e-10, rtol = 1e-10)
    sum(sol[-1, 2])
  }, c(log(ke0))), error = function(e) e)
say("RTMBode taped  : ", if (inherits(ad, "error"))
  paste("ERROR -", conditionMessage(ad)) else "ok")
say("=> a branch on time inside `dynamics` is REFUSED, not frozen: ",
    "func2tape() tapes the derivative function in both branches of ",
    "RTMBode::ode(), so `t` is an advector and RTMB blocks the ",
    "comparison. Loud, which is the good outcome.")

# smooth arithmetic in t is fine: it is only control flow that is barred
dyn_smooth <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(p[2] * exp(-p[3] * t) - p[1] * y[1]))
}
sm <- tryCatch(
  MakeTape(function(th) {
    "c" <- ADoverload("c")
    sum(RTMBode::ode(y = c(0), times = c(0, obs), func = dyn_smooth,
                     parms = c(exp(th[1]), R_inf, 0.5), method = "lsoda",
                     atol = 1e-10, rtol = 1e-10)[-1, 2])
  }, c(log(ke0))), error = function(e) e)
say("smooth f(t) in the derivative: ", if (inherits(sm, "error"))
  paste("ERROR -", conditionMessage(sm))
  else paste("ok, value", format(sm(c(log(ke0))), digits = 8)))

# --------------------------------------------------------------------
hr("2. approxfun() forcing inside dynamics")

ftab <- approxfun(c(0, Tinf - 1e-9, Tinf, 100), c(R_inf, R_inf, 0, 0),
                  rule = 2)
dyn_forcing <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(ftab(t) - p[1] * y[1]))
}
r2de <- tryCatch(
  deSolve::ode(y = c(0), times = c(0, obs), func = dyn_forcing,
               parms = c(ke0), method = "lsoda",
               atol = 1e-10, rtol = 1e-10)[-1, 2],
  error = function(e) e)
say("plain deSolve: ", if (inherits(r2de, "error"))
  paste("ERROR", conditionMessage(r2de))
  else paste("max err", format(max(abs(r2de - ana)), digits = 4)))
r2num <- tryCatch(
  RTMBode::ode(y = c(0), times = c(0, obs), func = dyn_forcing,
               parms = c(ke0), method = "lsoda",
               atol = 1e-10, rtol = 1e-10)[-1, 2],
  error = function(e) e)
say("RTMBode numeric: ", if (inherits(r2num, "error"))
  paste("ERROR", conditionMessage(r2num))
  else paste("max err", format(max(abs(r2num - ana)), digits = 4)))
r2ad <- tryCatch(
  MakeTape(function(th) {
    "c" <- ADoverload("c")
    sol <- RTMBode::ode(y = c(0), times = c(0, obs), func = dyn_forcing,
                        parms = c(exp(th[1])), method = "lsoda",
                        atol = 1e-10, rtol = 1e-10)
    sum(sol[-1, 2])
  }, c(log(ke0))), error = function(e) e)
if (inherits(r2ad, "error")) {
  say("taped: ERROR ", conditionMessage(r2ad))
} else {
  say("taped: value ", format(r2ad(c(log(ke0))), digits = 10),
      "  analytic ", format(sum(ana), digits = 10),
      "  frozen-at-0 ", format(sum(never_off), digits = 10))
}

# --------------------------------------------------------------------
hr("3. deSolve forcings= argument, through RTMBode")

dyn_plain <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}
r3 <- tryCatch(
  suppressWarnings(
    RTMBode::ode(y = c(100), times = c(0, obs), func = dyn_plain,
                 parms = c(ke0), method = "lsoda",
                 initforc = "forcc",
                 forcings = list(cbind(c(0, Tinf, 100), c(R_inf, 0, 0))))),
  error = function(e) e)
say(if (inherits(r3, "error"))
  paste("forcings= rejected:", conditionMessage(r3))
  else paste("forcings= did not error, but RTMBode's compiled shim",
             "'desolve_derivs' has no forcing hook, so the forcing is",
             "never read. First state:", format(r3[2, 2], digits = 8),
             "vs plain decay", format(100 * exp(-ke0 * obs[1]),
                                      digits = 8)))

# --------------------------------------------------------------------
hr("4. THE ROUTE THAT WORKS: rate as an extra parameter, per segment")

# The dynamics are wrapped so that the last nstate parameters are a
# constant rate vector added to the derivative. The rate is a genuine
# tape input, so it is differentiated exactly; it is held constant over
# one segment of a segmented solve, and the segment boundaries are data,
# so no branch on time is ever taken on the tape.
dyn_base <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}
wrap_rate <- function(dyn, np, ns) function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  r <- dyn(t, y, p[seq_len(np)])
  r <- if (is.list(r)) r[[1]] else r
  list(r + p[np + seq_len(ns)])
}
dyn_r <- wrap_rate(dyn_base, 1L, 1L)

seg_infusion <- function(th, brk, rates, want) {
  "c" <- ADoverload("c"); "[<-" <- ADoverload("[<-")
  ke <- exp(th[1])
  out <- numeric(length(want))
  y <- c(0)
  for (s in seq_len(length(brk) - 1L)) {
    a <- brk[s]; b <- brk[s + 1L]
    ts <- want[want > a & want <= b]
    gg <- c(a, ts, if (!length(ts) || ts[length(ts)] < b) b)
    sol <- RTMBode::ode(y = y, times = gg, func = dyn_r,
                        parms = c(ke, rates[s]), method = "lsoda",
                        atol = 1e-10, rtol = 1e-10)
    if (length(ts)) out[match(ts, want)] <- sol[seq_along(ts) + 1L, 2]
    y <- c(sol[nrow(sol), 2])
  }
  out
}

brk <- c(0, Tinf, max(obs))
tp4 <- MakeTape(function(th)
  seg_infusion(th, brk, c(R_inf, 0), obs), c(log(ke0)))
v4 <- as.numeric(tp4(c(log(ke0))))
say("max |segmented - analytic| = ", format(max(abs(v4 - ana)), digits = 4))
tp4s <- MakeTape(function(th)
  sum(seg_infusion(th, brk, c(R_inf, 0), obs)), c(log(ke0)))
report("d/d(log ke), infusion via rate parameter",
       as.numeric(tp4s$jacfun()(c(log(ke0)))),
       cfd(function(th) sum(inf_analytic(obs, R_inf, exp(th[1]), Tinf)),
           c(log(ke0))))

# --------------------------------------------------------------------
hr("5. infusion rate itself estimated (rate = amount / duration)")

tp5 <- MakeTape(function(th)
  sum(seg_infusion(th[1], brk, c(exp(th[2]) / Tinf, 0), obs)),
  c(log(ke0), log(R_inf * Tinf)))
report("d/d(log ke, log amount)",
       as.numeric(tp5$jacfun()(c(log(ke0), log(R_inf * Tinf)))),
       cfd(function(th) sum(inf_analytic(obs, exp(th[2]) / Tinf,
                                         exp(th[1]), Tinf)),
           c(log(ke0), log(R_inf * Tinf))))

# --------------------------------------------------------------------
hr("6. cost: one solve vs one solve per segment")

nrep <- 200
tp_one <- MakeTape(function(th) {
  "c" <- ADoverload("c")
  sum(RTMBode::ode(y = c(100), times = c(0, obs), func = dyn_base,
                   parms = c(exp(th[1])), method = "lsoda",
                   atol = 1e-8, rtol = 1e-8)[-1, 2])
}, c(log(ke0)))
j_one <- tp_one$jacfun(); j_seg <- tp4s$jacfun()
# ODEadjoint caches the solution keyed on the input vector, so the point
# must move or the benchmark measures the cache
xs <- log(ke0) + seq(-0.05, 0.05, length.out = nrep)
invisible(j_one(xs[1])); invisible(j_seg(xs[1]))
t_one <- system.time(for (i in seq_len(nrep))
  j_one(xs[i]))[["elapsed"]]
t_seg <- system.time(for (i in seq_len(nrep))
  j_seg(xs[i]))[["elapsed"]]
say("1 segment,  ", nrep, " gradients: ", format(t_one, digits = 3),
    "s (", format(1000 * t_one / nrep, digits = 3), " ms each)")
say("2 segments, ", nrep, " gradients: ", format(t_seg, digits = 3),
    "s (", format(1000 * t_seg / nrep, digits = 3), " ms each)")
say("ratio = ", format(t_seg / t_one, digits = 3), "x for 2x the solves")

# and the tape-build cost, which is paid once per fit
t_build <- system.time(for (i in seq_len(20))
  MakeTape(function(th) sum(seg_infusion(th, brk, c(R_inf, 0), obs)),
           c(log(ke0))))[["elapsed"]]
say("tape build, 20 x 2-segment: ", format(t_build, digits = 3), "s")

say("\ndone.")
