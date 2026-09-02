# Probe H1: does RTMBode::ode() pass deSolve `events` through, and does
# the numeric trajectory match the analytic multi-dose superposition?
#
# Model: one-compartment oral, repeated bolus doses into the depot.
#   dA/dt = -ka A ; dC/dt = ka A / V - ke C
# A dose D at time td adds D to A. By linearity the concentration is the
# superposition of single-dose curves shifted to each dose time:
#   C(t) = sum_{td <= t} D ka / (V (ka - ke)) (exp(-ke u) - exp(-ka u)),
#          u = t - td
# That closed form is an EXACT reference, so any mismatch is the solver
# or the event machinery, not a modelling choice.
#
# Also pins down the deSolve event-time contract (task item c):
#   - must an event time appear in `times`?
#   - is the reported value at a dose time pre- or post-dose?
#   - duplicated times at a dose instant.

suppressPackageStartupMessages({
  library(RTMB)
  library(RTMBode)
})
source("dev/ode/pk-common.R")

say <- function(...) cat(..., "\n", sep = "")
hr <- function(t) say("\n===== ", t, " =====")

ka <- 1.0; ke <- 0.2; V <- 10; D <- 100
dose_times <- c(0, 6, 12, 18)

multi_dose_analytic <- function(t, ka, ke, V, D, dose_times) {
  vapply(t, function(ti) {
    u <- ti - dose_times[dose_times <= ti]
    if (!length(u)) return(0)
    sum(D * ka / (V * (ka - ke)) * (exp(-ke * u) - exp(-ka * u)))
  }, 0)
}

# --------------------------------------------------------------------
hr("1. events= reaches deSolve at all (numeric path)")

obs <- c(1, 3, 5, 7, 9, 11, 13, 17, 20, 24)
# Dose at t=0 is the initial condition; later doses are events.
ev_times <- dose_times[-1]
grid <- sort(unique(c(0, obs, ev_times)))
ev <- data.frame(var = 1, time = ev_times, value = D, method = "add")

res <- tryCatch(
  RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn,
               parms = c(ka, ke, V), method = "lsoda",
               atol = 1e-10, rtol = 1e-10,
               events = list(data = ev)),
  error = function(e) e)
if (inherits(res, "error")) {
  say("FAIL events rejected: ", conditionMessage(res))
} else {
  say("OK   solve returned ", nrow(res), " x ", ncol(res))
  num <- res[match(obs, res[, 1]), 3]
  ana <- multi_dose_analytic(obs, ka, ke, V, D, dose_times)
  say("max |ode - analytic| = ", format(max(abs(num - ana)), digits = 4))
  say("max rel             = ",
      format(max(abs(num - ana) / pmax(abs(ana), 1e-8)), digits = 4))
  print(data.frame(t = obs, ode = num, analytic = ana))
}

# --------------------------------------------------------------------
hr("2. control: same solve with NO events (single dose)")

res0 <- RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn,
                     parms = c(ka, ke, V), method = "lsoda",
                     atol = 1e-10, rtol = 1e-10)
num0 <- res0[match(obs, res0[, 1]), 3]
ana0 <- multi_dose_analytic(obs, ka, ke, V, D, 0)
say("max |ode - analytic| single dose = ",
    format(max(abs(num0 - ana0)), digits = 4))
say("events changed the trajectory: ",
    if (!inherits(res, "error")) max(abs(num - num0)) > 1 else NA)

# --------------------------------------------------------------------
hr("3. contract: must event times be in `times`?")

grid_no_ev <- sort(unique(c(0, obs)))   # 6, 12, 18 absent
r3 <- tryCatch(
  RTMBode::ode(y = c(D, 0), times = grid_no_ev, func = pk_dyn,
               parms = c(ka, ke, V), method = "lsoda",
               atol = 1e-10, rtol = 1e-10,
               events = list(data = ev)),
  error = function(e) e, warning = function(w) w)
if (inherits(r3, "condition")) {
  say(class(r3)[1], ": ", conditionMessage(r3))
} else {
  n3 <- r3[match(obs, r3[, 1]), 3]
  a3 <- multi_dose_analytic(obs, ka, ke, V, D, dose_times)
  say("solved without event times in grid; max err = ",
      format(max(abs(n3 - a3)), digits = 4))
}

# --------------------------------------------------------------------
hr("4. contract: observation exactly AT a dose time - pre or post?")

obs4 <- c(5.999999, 6, 6.000001)
grid4 <- sort(unique(c(0, obs4, ev_times)))
r4 <- RTMBode::ode(y = c(D, 0), times = grid4, func = pk_dyn,
                   parms = c(ka, ke, V), method = "lsoda",
                   atol = 1e-10, rtol = 1e-10,
                   events = list(data = ev))
show <- r4[r4[, 1] %in% obs4, , drop = FALSE]
print(show)
say("depot just before 6: ", format(show[1, 2], digits = 8))
say("depot at         6: ", format(show[2, 2], digits = 8))
say("depot just after  6: ", format(show[3, 2], digits = 8))
say("=> value AT the dose time is ",
    if (abs(show[2, 2] - show[1, 2]) > 1) "POST-dose" else "PRE-dose")

# --------------------------------------------------------------------
hr("5. contract: duplicated time at the dose instant")

grid5 <- c(0, 3, 6, 6, 9)
r5 <- tryCatch(
  RTMBode::ode(y = c(D, 0), times = grid5, func = pk_dyn,
               parms = c(ka, ke, V), method = "lsoda",
               atol = 1e-10, rtol = 1e-10,
               events = list(data = ev)),
  error = function(e) e, warning = function(w) w)
if (inherits(r5, "condition")) {
  say(class(r5)[1], ": ", conditionMessage(r5))
} else {
  print(r5)
  say("both rows at t=6 identical: ",
      isTRUE(all.equal(r5[3, ], r5[4, ])))
}

# --------------------------------------------------------------------
hr("6. contract: event var by NAME vs by index")

y_named <- c(depot = D, central = 0)
ev_name <- data.frame(var = "depot", time = ev_times, value = D,
                      method = "add")
r6 <- tryCatch(
  RTMBode::ode(y = y_named, times = grid, func = pk_dyn,
               parms = c(ka, ke, V), method = "lsoda",
               atol = 1e-10, rtol = 1e-10,
               events = list(data = ev_name)),
  error = function(e) e)
if (inherits(r6, "error")) {
  say("FAIL name lookup: ", conditionMessage(r6))
} else {
  say("name lookup OK; matches index form: ",
      isTRUE(all.equal(unname(r6[, 3]), unname(res[, 3]))))
}

# --------------------------------------------------------------------
hr("7. methods replace / multiply, numerically")

ev_rep <- data.frame(var = 1, time = 6, value = 50, method = "replace")
r7 <- RTMBode::ode(y = c(D, 0), times = c(0, 5.999, 6, 6.001, 8),
                   func = pk_dyn, parms = c(ka, ke, V), method = "lsoda",
                   atol = 1e-10, rtol = 1e-10,
                   events = list(data = ev_rep))
print(r7)
ev_mul <- data.frame(var = 1, time = 6, value = 2, method = "multiply")
r8 <- RTMBode::ode(y = c(D, 0), times = c(0, 5.999, 6, 6.001, 8),
                   func = pk_dyn, parms = c(ka, ke, V), method = "lsoda",
                   atol = 1e-10, rtol = 1e-10,
                   events = list(data = ev_mul))
print(r8)

# --------------------------------------------------------------------
hr("8. integrator sensitivity: does every adaptive method honour events?")

meths <- c("lsoda", "lsode", "vode", "radau", "bdf", "adams",
           "impAdams", "ode45", "ode23", "daspk", "lsodes", "lsodar")
for (m in meths) {
  rr <- tryCatch(
    suppressWarnings(RTMBode::ode(y = c(D, 0), times = grid, func = pk_dyn,
                 parms = c(ka, ke, V), method = m,
                 atol = 1e-10, rtol = 1e-10,
                 events = list(data = ev))),
    error = function(e) e)
  if (inherits(rr, "error")) {
    say(sprintf("%-10s ERROR %s", m, conditionMessage(rr)))
  } else if (nrow(rr) != length(grid)) {
    say(sprintf("%-10s SHORT %d of %d rows", m, nrow(rr), length(grid)))
  } else {
    nn <- rr[match(obs, rr[, 1]), 3]
    aa <- multi_dose_analytic(obs, ka, ke, V, D, dose_times)
    say(sprintf("%-10s max err %.3e", m, max(abs(nn - aa))))
  }
}

say("\ndone.")
