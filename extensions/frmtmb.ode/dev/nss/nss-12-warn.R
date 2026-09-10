# The false-alarm rate of the run-in warning, on the schedule space the
# field actually writes.
#
# A warning that fires on a model whose answer is right is a real cost,
# so the classification is against the TRUTH and not against itself:
# frm_lincmt() at its default n_ss = Inf is the exact steady state, so
# the error of what frm_ode() returns over one dosing interval is known
# for every row. A row is
#
#   correct silence  no warning, error <= ss_tol
#   correct warning  warning,    error >  ss_tol
#   FALSE ALARM      warning,    error <= ss_tol
#   MISS             no warning, error >  ss_tol
#
# Run with NSS_ARM=ref for the base commit, where `ss_extrapolate` is
# absorbed by `...` and the warning reports the cycle-to-cycle movement.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-12-warn.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()
ARM <- Sys.getenv("NSS_ARM", "lane")
SS_TOL <- 1e-6

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
one_iv <- function(t, y, p) list(c(-p[1L] * y[1L]))
two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}
two_iv <- function(t, y, p) {
  list(c(-(p[1L] + p[2L]) * y[1L] + p[3L] * y[2L],
         p[2L] * y[1L] - p[3L] * y[2L]))
}

cases <- list()
add <- function(drug, ncmt, depot, ke, k12, k21, ka, ii, dur = 0,
                nss = 20L, tol = 1e-8) {
  cases[[length(cases) + 1L]] <<- list(
    drug = drug, ncmt = ncmt, depot = depot, ke = ke, k12 = k12,
    k21 = k21, ka = ka, ii = ii, dur = dur, nss = nss, tol = tol)
}
# half-lives and intervals from ordinary therapeutics: a short-acting
# oral analgesic, a twice-daily antibiotic, a daily statin, a weekly
# bisphosphonate, a monoclonal antibody dosed every three weeks
add("1 cmt oral, t1/2 3.5 h, q8h",     1, TRUE,  0.2,   NA, NA, 1.1, 8)
add("1 cmt oral, t1/2 3.5 h, q12h",    1, TRUE,  0.2,   NA, NA, 1.1, 12)
add("1 cmt oral, t1/2 14 h, q24h",     1, TRUE,  0.05,  NA, NA, 1.1, 24)
add("1 cmt oral, t1/2 35 h, q24h",     1, TRUE,  0.02,  NA, NA, 1.0, 24)
add("1 cmt oral, t1/2 69 h, q12h",     1, TRUE,  0.01,  NA, NA, 1.0, 12)
add("1 cmt oral, t1/2 139 h, q12h",    1, TRUE,  0.005, NA, NA, 1.0, 12)
add("1 cmt iv,   t1/2 3.5 h, q8h",     1, FALSE, 0.2,   NA, NA, NA,  8)
add("1 cmt iv,   t1/2 35 h, q24h",     1, FALSE, 0.02,  NA, NA, NA,  24)
add("2 cmt oral, t1/2z 23 h, q8h",     2, TRUE,  0.2,  0.4,  0.1,   1.1, 8)
add("2 cmt oral, t1/2z 23 h, q24h",    2, TRUE,  0.2,  0.4,  0.1,   1.1, 24)
add("2 cmt oral, t1/2z 42 h, q12h",    2, TRUE,  0.5,  1.0,  0.05,  1.1, 12)
add("2 cmt oral, t1/2z 107 h, q24h",   2, TRUE,  0.15, 0.3,  0.02,  1.0, 24)
add("2 cmt oral, t1/2z 107 h, q12h",   2, TRUE,  0.15, 0.3,  0.02,  1.0, 12)
add("2 cmt oral, t1/2z 265 h, q24h",   2, TRUE,  0.1,  0.2,  0.008, 1.0, 24)
add("2 cmt oral, t1/2z 670 h, q24h",   2, TRUE,  0.08, 0.15, 0.003, 1.0, 24)
add("2 cmt oral, t1/2z 2089 h, q168h", 2, TRUE,  0.05, 0.1,  0.001, 1.0, 168)
add("2 cmt iv,   t1/2z 23 h, q8h",     2, FALSE, 0.2,  0.4,  0.1,   NA,  8)
add("2 cmt iv,   t1/2z 107 h, q24h",   2, FALSE, 0.15, 0.3,  0.02,  NA,  24)
add("2 cmt iv,   t1/2z 265 h, q24h",   2, FALSE, 0.1,  0.2,  0.008, NA,  24)
# infusions: a one-hour iv and a slow-release oral read as zero order
add("2 cmt iv 1 h infusion, t1/2z 23 h, q8h",   2, FALSE, 0.2, 0.4, 0.1,  NA, 8,  1)
add("2 cmt iv 2 h infusion, t1/2z 107 h, q24h", 2, FALSE, 0.15, 0.3, 0.02, NA, 24, 2)
add("1 cmt iv 4 h infusion, t1/2 35 h, q24h",   1, FALSE, 0.02, NA, NA,   NA, 24, 4)
# FLIP-FLOP: ka near lambda_z, which is what an extended-release or
# depot formulation is written to produce. Every row above runs ka at
# 1.0 or 1.1 per hour, so none of them can enter this regime, and it is
# the regime where the correction can LOSE to truncation
# (dev/reviews/2026-09-10-nss.md NIT 3).
add("2 cmt oral FLIP-FLOP, t1/2z 333 h, q24h", 2, TRUE, 0.0223, 0.0445, 0.0067, 0.0031, 24)
add("1 cmt oral flip-flop, t1/2 69 h, q24h",   1, TRUE, 0.01,   NA,     NA,     0.008,  24)
add("1 cmt oral flip-flop, t1/2 17 h, q12h",   1, TRUE, 0.04,   NA,     NA,     0.03,   12)
# a run-in far too short, where `keep` is not even full
add("2 cmt oral 107 h q24, n_ss = 1",  2, TRUE, 0.15, 0.3, 0.02, 1.0, 24, 0, 1L)
add("2 cmt oral 107 h q24, n_ss = 2",  2, TRUE, 0.15, 0.3, 0.02, 1.0, 24, 0, 2L)
add("2 cmt oral 107 h q24, n_ss = 4",  2, TRUE, 0.15, 0.3, 0.02, 1.0, 24, 0, 4L)
# LARGE n_ss in the stand-down band, which `?frm_ode` sends a
# long-half-life user to and where the warning used to go silent
# (dev/reviews/2026-09-10-nss.md NIT 11). r = 0.98999, above rcap.
# These carry atol = rtol = 1e-12 rather than the shipped 1e-8. A
# run-in 1000 cycles long chains 1001 solves, so at 1e-8 the
# integrator contributes about 2e-05 and the run-in's own shortfall is
# no longer the error being classified. That is a real limit on
# "raise n_ss" and `?frm_ode` now says so.
add("1 cmt iv t1/2 1653 h q24, n_ss = 650",  1, FALSE, log(2)/1653, NA, NA, NA, 24, 0, 650L,  1e-12)
add("1 cmt iv t1/2 1653 h q24, n_ss = 1000", 1, FALSE, log(2)/1653, NA, NA, NA, 24, 0, 1000L, 1e-12)
add("1 cmt iv t1/2 1653 h q24, n_ss = 2000", 1, FALSE, log(2)/1653, NA, NA, NA, 24, 0, 2000L, 1e-12)

exact_of <- function(cs, tt) {
  ev <- data.frame(time = 0,
                   state = if (cs$depot) "depot" else "central",
                   value = 100, ii = cs$ii, addl = 0L, ss = TRUE)
  if (cs$dur > 0) ev$duration <- cs$dur
  P <- if (cs$ncmt == 1L)
    list(ke = cs$ke, ka = cs$ka, V = 10)
  else list(ke = cs$ke, k12 = cs$k12, k21 = cs$k21, ka = cs$ka, V = 10)
  P <- P[!vapply(P, function(z) is.na(z[1L]), TRUE)]
  as.numeric(frm_lincmt(parms = P, times = tt, ncmt = cs$ncmt,
                        depot = cs$depot, events = ev))
}
ode_of <- function(cs, tt, ext) {
  dyn <- if (cs$ncmt == 1L && cs$depot) one_oral
         else if (cs$ncmt == 1L) one_iv
         else if (cs$depot) two_oral else two_iv
  ns <- if (cs$ncmt == 1L) (if (cs$depot) 2L else 1L)
        else (if (cs$depot) 3L else 2L)
  out <- if (cs$depot) 2L else 1L
  pv <- if (cs$ncmt == 1L) list(cs$ke, cs$ka) else
    list(cs$ke, cs$k12, cs$k21, cs$ka)
  pv <- pv[!vapply(pv, function(z) is.na(z[1L]), TRUE)]
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = cs$ii,
                   ss = TRUE)
  if (cs$dur > 0) ev$duration <- cs$dur
  args <- list(dyn, init = rep(list(0), ns), times = tt, parms = pv,
               events = ev, output = out, n_ss = cs$nss,
               ss_tol = SS_TOL, atol = cs$tol, rtol = cs$tol)
  if (ARM != "ref") args$ss_extrapolate <- ext
  warned <- FALSE
  v <- withCallingHandlers(do.call(frm_ode, args),
                           warning = function(w) {
                             if (grepl("run-in", conditionMessage(w)))
                               warned <<- TRUE
                             invokeRestart("muffleWarning")
                           })
  list(v = as.numeric(v) / 10, warned = warned)
}

arms <- if (ARM == "ref") list(base = NA) else
  list(extrapolated = TRUE, truncated = FALSE)
for (an in names(arms)) {
  cat(sprintf("\n=== %s (%s) ===\n", an, ARM))
  cat(sprintf("%-42s %11s %9s %s\n", "case", "error", "warned",
              "verdict"))
  tally <- c(silence = 0L, warn = 0L, false = 0L, miss = 0L)
  for (cs in cases) {
    tt <- seq(0, cs$ii, length.out = 25)
    ref <- exact_of(cs, tt)
    z <- ode_of(cs, tt, arms[[an]])
    err <- max(abs(z$v - ref)) / max(abs(ref))
    bad <- err > SS_TOL
    verdict <- if (z$warned && bad) "correct warning"
               else if (!z$warned && !bad) "correct silence"
               else if (z$warned) "FALSE ALARM" else "MISS"
    k <- c("correct silence" = "silence", "correct warning" = "warn",
           "FALSE ALARM" = "false", "MISS" = "miss")[[verdict]]
    tally[[k]] <- tally[[k]] + 1L
    cat(sprintf("%-42s %11.3e %9s %s\n", cs$drug, err, z$warned,
                verdict))
  }
  cat(sprintf("  %d cases: %d correct silence, %d correct warning, ",
              length(cases), tally[["silence"]], tally[["warn"]]))
  cat(sprintf("%d FALSE ALARM, %d MISS\n", tally[["false"]],
              tally[["miss"]]))
}
