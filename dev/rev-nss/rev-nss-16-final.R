# REVIEW of lane nss round 2. Three scoped items.
#
# 1. THE CURVATURE CONTROL. The lane says the residual gap at `rlow`
#    and `rcap` is the probe's own `f''(r) * eps` and not a corner,
#    checked against `2 / (1 - r)^3`. A third junction is not the right
#    third test, because the other two junctions are not of this shape.
#    The right one is a NON-junction: if the probe reads curvature, it
#    must obey the same law where there is provably nothing to find,
#    and `F(r)` is exactly `r / (1 - r)` on the whole open interval
#    between the gate and the cap. So the law is checked at five
#    interior points as well as at the two junctions, over five
#    brackets rather than four.
#
# 2. THE WARNING, from scratch. `undone` is new arithmetic on a
#    diagnostic that has failed open twice, so the absent cases are
#    rebuilt rather than re-run:
#      2a the declined part nonzero where the true error is not, which
#         would be a false alarm. The gate band is where `fac` is below
#         `whole` by construction, and `ode_ss_trust` = 0.5 splits it
#         into a modelled half and a movement-fallback half;
#      2b the correction applied IN FULL, so `undone` is exactly zero,
#         and the answer still wrong. That is the flip-flop and
#         two-near-modes family, where `r` is an honest-looking
#         contraction and the geometric model is not the truth;
#      2c the boundary where it goes quiet, between the lane's 1000 and
#         2000 rows.
#
# Script path: dev/rev-nss/rev-nss-16-final.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()
ext <- frmtmb.ode:::ode_ss_extrapolate
RCAP <- frmtmb.ode:::ode_ss_rcap
RLOW <- frmtmb.ode:::ode_ss_rlow
TRUST <- frmtmb.ode:::ode_ss_trust
SS_TOL <- 1e-6
cat("\nrcap", RCAP, " rlow", RLOW, " trust", TRUST, "\n")

## ---- 1. the curvature control, at non-junctions too ---------------
mk <- function(r) c(0, 1 / r, 1 / r + 1)
FF <- function(r) {
  if (r == 0) return(0)
  v <- mk(r)
  as.numeric(ext(v[1L], v[2L], v[3L], 1e-14, 1e-14)$y - v[3L])
}
one_sided_gap <- function(x0, eps) {
  l <- (FF(x0) - FF(x0 - eps)) / eps
  r <- (FF(x0 + eps) - FF(x0)) / eps
  abs(l - r)
}
cat("\n== 1. gap / eps against 2 / (1 - r)^3 ==\n")
cat("   the two junctions, and five points where there is no junction",
    "at all\n\n")
pts <- c(RLOW, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, RCAP)
cat(sprintf("%9s %9s", "r", "2/(1-r)^3"))
eps <- 10^-(2:6)
for (e in eps) cat(sprintf(" %12s", format(e)))
cat("\n")
for (x0 in pts) {
  cat(sprintf("%9.4f %9.4g", x0, 2 / (1 - x0)^3))
  for (e in eps) cat(sprintf(" %12.5g", one_sided_gap(x0, e) / e))
  cat(if (x0 %in% c(RLOW, RCAP)) "   <- junction\n" else "\n")
}
cat("\n   ratio of the eps = 1e-5 column to 2/(1-r)^3:\n")
for (x0 in pts)
  cat(sprintf("     r = %-7.4f  %8.5f\n", x0,
              (one_sided_gap(x0, 1e-5) / 1e-5) / (2 / (1 - x0)^3)))
cat(sprintf("\n   gate overshoot: max F / (r/(1-r)) on (0, rlow] = %.4f\n",
            max(vapply(seq(1e-4, RLOW, length.out = 2000),
                       function(r) FF(r) / (r / (1 - r)), 0))))

## ---- 2. the warning, absent cases rebuilt --------------------------
one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
go <- function(ncmt, pv, ii, n_ss, e, tol) {
  dyn <- if (ncmt == 1L) one_oral else two_oral
  ns <- if (ncmt == 1L) 2L else 3L
  tt <- seq(0, ii, length.out = 25)
  w <- character(0)
  v <- withCallingHandlers(
    frm_ode(dyn, init = rep(list(0), ns), times = tt, parms = pv,
            events = data.frame(time = 0, state = 1L, value = 100,
                                ii = ii, ss = TRUE),
            output = 2L, n_ss = n_ss, ss_tol = SS_TOL,
            ss_extrapolate = e, atol = tol, rtol = tol),
    warning = function(z) { w <<- c(w, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  ww <- w[grepl("run-in", w)]
  list(v = as.numeric(v) / 10, warned = length(ww) > 0L,
       said = if (!length(ww)) NA_real_ else
         as.numeric(sub("^.*is still ([0-9.e+-]+) [(]relative.*$",
                        "\\1", ww[[1L]])))
}
lin <- function(ncmt, P, ii)
  as.numeric(frm_lincmt(parms = P, times = seq(0, ii, length.out = 25),
                        ncmt = ncmt, depot = TRUE,
                        events = data.frame(time = 0, state = "depot",
                                            value = 100, ii = ii,
                                            addl = 0L, ss = TRUE)))
rows <- list()
say <- function(tag, ncmt, pv, P, ii, n_ss, tol = 1e-10) {
  ref <- lin(ncmt, P, ii)
  for (e in c(TRUE, FALSE)) {
    z <- go(ncmt, pv, ii, n_ss, e, tol)
    err <- max(abs(z$v - ref)) / max(abs(ref))
    bad <- err > SS_TOL
    verdict <- if (z$warned && bad) "ok warn"
               else if (!z$warned && !bad) "ok silent"
               else if (z$warned) "FALSE ALARM" else "MISS"
    rows[[length(rows) + 1L]] <<- data.frame(
      tag = tag, n_ss = n_ss, arm = if (e) "TRUE" else "FALSE",
      tol = tol, err = err, said = z$said, verdict = verdict)
  }
}

# 2a. the gate band, where `fac` is below `whole` by construction.
# r = exp(-ke * ii) is set directly; short run-ins are the only place
# the differences survive the damping at these ratios.
cat("\n== 2a. the gate band: is `undone` a false alarm? ==\n")
for (rt in c(0.010, 0.020, 0.035, 0.049, 0.06, 0.12)) {
  ke <- -log(rt) / 24
  for (n in c(3L, 4L, 5L, 6L, 8L))
    say(sprintf("gate band r=%.3f", rt), 1L, list(ke, 1.0),
        list(ke = ke, ka = 1.0, V = 10), 24, n)
}
# 2b. the correction applied IN FULL, so `undone` is exactly 0, and the
# answer still wrong: flip-flop and two near-equal modes.
cat("== 2b. undone == 0 but the answer wrong ==\n")
say("flip-flop ka/lz 1.5", 2L, list(0.022267, 0.044534, 0.006680,
                                    0.003125),
    list(ke = 0.022267, k12 = 0.044534, k21 = 0.006680,
         ka = 0.003125, V = 10), 24, 20L)
say("flip-flop ka/lz 0.8", 2L, list(0.089068, 0.178135, 0.026720,
                                    0.006667),
    list(ke = 0.089068, k12 = 0.178135, k21 = 0.026720,
         ka = 0.006667, V = 10), 24, 20L)
say("two near modes", 2L, list(0.10, 0.02, 0.0805, 1.0),
    list(ke = 0.10, k12 = 0.02, k21 = 0.0805, ka = 1.0, V = 10),
    24, 20L)
say("big peripheral 100x", 2L, list(0.05, 0.50, 0.005, 1.0),
    list(ke = 0.05, k12 = 0.5, k21 = 0.005, ka = 1, V = 10), 24, 20L)
say("big peripheral 100x", 2L, list(0.05, 0.50, 0.005, 1.0),
    list(ke = 0.05, k12 = 0.5, k21 = 0.005, ka = 1, V = 10), 24, 320L)
# 2c. the boundary where it goes quiet
cat("== 2c. the 1000 to 2000 boundary ==\n")
KE <- log(2) / 1653
for (n in c(650L, 1000L, 1200L, 1400L, 1700L, 2000L))
  say("t1/2 1653 h", 1L, list(KE, 1.0),
      list(ke = KE, ka = 1.0, V = 10), 24, n, 1e-12)

d <- do.call(rbind, rows)
cat(sprintf("\n%-24s %6s %6s %7s %11s %11s %10s %s\n", "case", "n_ss",
            "arm", "tol", "true error", "it said", "said/true",
            "verdict"))
for (i in seq_len(nrow(d)))
  cat(sprintf("%-24s %6d %6s %7.0e %11.3e %11s %10s %s\n", d$tag[i],
              d$n_ss[i], d$arm[i], d$tol[i], d$err[i],
              if (is.na(d$said[i])) "-" else
                format(signif(d$said[i], 3)),
              if (is.na(d$said[i]) || d$err[i] == 0) "-" else
                format(signif(d$said[i] / d$err[i], 3)),
              d$verdict[i]))
cat("\n-- tally over", nrow(d), "runs --\n")
print(table(d$arm, d$verdict))
cat("\ndone\n")
