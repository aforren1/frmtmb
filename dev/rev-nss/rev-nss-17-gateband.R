# REVIEW of lane nss round 2: are the gate-band false alarms new?
#
# `dev/rev-nss/rev-nss-16.log` finds 13 false alarms in 82 runs, all
# with the measured ratio inside the gate, [0, ode_ss_rlow], and all at
# short run-ins. They are not misses and they are not wrong answers,
# but the lane claims 0 false alarms, so the question is whether the
# `undone` term of round 2 put them there.
#
# The mechanism says it did. In the gate band `fac` is `whole` times
# `gate`, and `gate` falls below `ode_ss_trust` = 0.5 once
# `r / rlow` is under about 0.72, so `modelled` is FALSE and both
# `full` and `undone` fall back to `move`. `move` OVERSTATES the
# distance to the limit by about 1 / r there, which is 100x at
# r = 0.01. That is the mirror of the defect round 0 removed: the
# movement understates when r is near 1 and overstates when r is near
# 0, and the fallback now quotes it at the near-0 end.
#
# Round 1 cannot do this: it had no `modelled` test, `bad` needed `r`
# OUTSIDE (0, 1), and r in the gate band is inside it. The `noundone`
# variant of the lane's own see-it-fail harness restores exactly that
# pair of lines, so it is the right control and it doubles as the
# second variant this review was asked to spot-check end to end.
#
# REV_ARM=noundone selects it.
#
# Script path: dev/rev-nss/rev-nss-17-gateband.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()
SS_TOL <- 1e-6

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
go <- function(ke, n_ss, e) {
  tt <- seq(0, 24, length.out = 25)
  w <- character(0)
  v <- withCallingHandlers(
    frm_ode(one_oral, init = list(0, 0), times = tt,
            parms = list(ke, 1.0),
            events = data.frame(time = 0, state = 1L, value = 100,
                                ii = 24, ss = TRUE),
            output = 2L, n_ss = n_ss, ss_tol = SS_TOL,
            ss_extrapolate = e, atol = 1e-10, rtol = 1e-10),
    warning = function(z) { w <<- c(w, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  ww <- w[grepl("run-in", w)]
  list(v = as.numeric(v) / 10, warned = length(ww) > 0L,
       said = if (!length(ww)) NA_real_ else
         as.numeric(sub(paste0("^.*is (?:still|about) ([0-9.e+-]+) ",
                               "[(]relative.*$"), "\\1", ww[[1L]])))
}
cat(sprintf("\n%-8s %5s %6s %11s %11s %10s %s\n", "r", "n_ss", "arm",
            "true error", "it said", "said/true", "verdict"))
tal <- character(0)
for (rt in c(0.010, 0.020, 0.035, 0.049, 0.060, 0.120)) {
  ke <- -log(rt) / 24
  ref <- as.numeric(frm_lincmt(
    parms = list(ke = ke, ka = 1.0, V = 10),
    times = seq(0, 24, length.out = 25), ncmt = 1, depot = TRUE,
    events = data.frame(time = 0, state = "depot", value = 100,
                        ii = 24, addl = 0L, ss = TRUE)))
  for (n in c(3L, 4L, 5L, 6L, 8L)) for (e in c(TRUE, FALSE)) {
    z <- go(ke, n, e)
    err <- max(abs(z$v - ref)) / max(abs(ref))
    bad <- err > SS_TOL
    v <- if (z$warned && bad) "ok warn"
         else if (!z$warned && !bad) "ok silent"
         else if (z$warned) "FALSE ALARM" else "MISS"
    tal <- c(tal, paste(if (e) "TRUE" else "FALSE", v))
    cat(sprintf("%-8.3f %5d %6s %11.3e %11s %10s %s\n", rt, n,
                if (e) "TRUE" else "FALSE", err,
                if (is.na(z$said)) "-" else format(signif(z$said, 3)),
                if (is.na(z$said)) "-" else
                  format(signif(z$said / err, 3)), v))
  }
}
cat("\n-- tally over", length(tal), "runs --\n")
print(table(tal))
cat("\ndone\n")
