# REVIEW of lane nss round 1, attack 1: fresh absent-cases for the
# rewritten report.
#
# Round 1 took the round-0 constructions verbatim, so re-running them
# tests nothing new. This sweep was built from the report's own code
# rather than from the old table, by asking where each quantity it can
# form has a hole.
#
#   rel <- max(move) / scale
#   if (three) {
#     rel <- TRUE ? max(resid)/scale : max(tail + resid)/scale
#     moving <- move > ss_tol * scale
#     bad    <- moving & !(r > 0 & r < 1)
#     if (any(bad)) rel <- max(rel, max(move[bad]) / scale)
#     stalled <- any(bad & r >= 1)
#   }
#
# Three holes are visible in that, and each gets its own family here.
#
# 1. THE RAMP BAND. `bad` needs `r` OUTSIDE (0, 1), so a state with
#    `r` between `ode_ss_rcap` and 1 never reaches the `move` fallback.
#    There `tail` is the RAMPED correction, which the shape deliberately
#    makes small: `dev/rev-nss/rev-nss-12.log` measures F(0.999) at
#    1.175 where the geometric tail is 999, so `tail` understates the
#    distance by 850x by construction. Long terminal half-lives.
#
# 2. `move` IS THE OLD QUANTITY. `moving` gates the fallback on the
#    cycle-to-cycle movement, which understates the distance to the
#    limit by 1 / (1 - r). That is the exact defect round 0 rewrote the
#    warning to remove, and it is back as a gate rather than as a
#    report.
#
# 3. `scale` IS SHARED ACROSS STATES. `move`, `tail` and `resid` are
#    per state and every one of them is divided by `max(abs(y))` over
#    ALL states, so a state much smaller than the largest cannot raise
#    `rel` past `ss_tol` however wrong it is. A two-compartment model
#    with a large peripheral compartment does this on its own.
#
# Classification is against the TRUTH, `frm_lincmt()` at n_ss = Inf,
# over 25 points of one dosing interval, never against the run itself.
#
# Script path: dev/rev-nss/rev-nss-13-report.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()
SS_TOL <- 1e-6

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))

cases <- list()
add <- function(tag, ncmt, ke, k12, k21, ka, ii, n_ss) {
  cases[[length(cases) + 1L]] <<- list(
    tag = tag, ncmt = ncmt, ke = ke, k12 = k12, k21 = k21, ka = ka,
    ii = ii, n_ss = n_ss)
}
# family 1: the ramp band. r = exp(-lambda_z * ii) above ode_ss_rcap
for (n in c(20L, 40L)) {
  add("ramp band, t1/2z 1653 h", 1, log(2) / 1653, NA, NA, 1.0, 24, n)
  add("ramp band, t1/2z 3306 h", 1, log(2) / 3306, NA, NA, 1.0, 24, n)
  add("ramp band, t1/2z 16530 h", 1, log(2) / 16530, NA, NA, 1.0, 24, n)
  add("ramp band, 2 cmt t1/2z 2000 h", 2, 0.05, 0.1, 0.00042, 1.0, 24, n)
}
# family 2: a big peripheral compartment, so the observed state is small
# against the shared scale
for (n in c(20L, 40L)) {
  add("big peripheral 20x, t1/2z 107 h", 2, 0.15, 0.30, 0.0200, 1, 24, n)
  add("big peripheral 100x, t1/2z 340 h", 2, 0.05, 0.50, 0.0050, 1, 24, n)
  add("big peripheral 400x, t1/2z 900 h", 2, 0.02, 0.80, 0.0020, 1, 24, n)
  add("big peripheral 2000x, t1/2z 2800 h", 2, 0.01, 2.0, 0.0010, 1, 24, n)
}
# family 3: short run-ins, where `keep` is not full
for (n in 1:6) {
  add("short run-in, 1 cmt t1/2 35 h", 1, 0.02, NA, NA, 1.0, 24, n)
  add("short run-in, 2 cmt t1/2z 107 h", 2, 0.15, 0.3, 0.02, 1.0, 24, n)
}
# family 4: flip-flop, where the geometric model is wrong but r looks
# like an honest contraction
for (n in c(20L, 40L)) {
  add("flip-flop ka/lz 1.5, t1/2z 333 h", 2, 0.02227, 0.04453,
      0.006680, 0.003125, 24, n)
  add("flip-flop ka/lz 0.8, t1/2z 83 h", 2, 0.08907, 0.17814,
      0.026720, 0.006667, 24, n)
}

run1 <- function(cs, ext) {
  dyn <- if (cs$ncmt == 1L) one_oral else two_oral
  ns <- if (cs$ncmt == 1L) 2L else 3L
  pv <- if (cs$ncmt == 1L) list(cs$ke, cs$ka) else
    list(cs$ke, cs$k12, cs$k21, cs$ka)
  ev <- data.frame(time = 0, state = 1L, value = 100, ii = cs$ii,
                   ss = TRUE)
  tt <- seq(0, cs$ii, length.out = 25)
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(frm_ode(dyn, init = rep(list(0), ns), times = tt,
                     parms = pv, events = ev, output = 2L,
                     n_ss = cs$n_ss, ss_tol = SS_TOL,
                     ss_extrapolate = ext, atol = 1e-10, rtol = 1e-10),
             error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(e) { w <<- c(w, conditionMessage(e))
                            invokeRestart("muffleWarning") })
  ww <- w[grepl("run-in", w)]
  list(v = if (is.character(v)) NA_real_ else as.numeric(v) / 10,
       warned = length(ww) > 0L,
       said = if (!length(ww)) NA_real_ else
         as.numeric(sub("^.*is still ([0-9.e+-]+) [(]relative.*$",
                        "\\1", ww[[1L]])))
}
truth <- function(cs) {
  P <- if (cs$ncmt == 1L) list(ke = cs$ke, ka = cs$ka, V = 10) else
    list(ke = cs$ke, k12 = cs$k12, k21 = cs$k21, ka = cs$ka, V = 10)
  as.numeric(frm_lincmt(parms = P,
                        times = seq(0, cs$ii, length.out = 25),
                        ncmt = cs$ncmt, depot = TRUE,
                        events = data.frame(time = 0, state = "depot",
                                            value = 100, ii = cs$ii,
                                            addl = 0L, ss = TRUE)))
}

cat(sprintf("\n%-34s %4s %5s %11s %8s %11s %9s %s\n", "case", "n_ss",
            "arm", "true error", "warned", "it said", "said/true",
            "verdict"))
tally <- list()
for (cs in cases) {
  ref <- truth(cs)
  for (ext in c(TRUE, FALSE)) {
    z <- run1(cs, ext)
    err <- if (all(is.na(z$v))) NA_real_ else
      max(abs(z$v - ref)) / max(abs(ref))
    bad <- is.finite(err) && err > SS_TOL
    verdict <- if (z$warned && bad) "ok warn"
               else if (!z$warned && !bad) "ok silent"
               else if (z$warned) "FALSE ALARM" else "MISS"
    k <- if (ext) "TRUE" else "FALSE"
    tally[[k]] <- c(tally[[k]], verdict)
    cat(sprintf("%-34s %4d %5s %11.3e %8s %11s %9s %s\n", cs$tag,
                cs$n_ss, k, err, z$warned,
                if (is.na(z$said)) "-" else format(signif(z$said, 3)),
                if (is.na(z$said) || !is.finite(err) || err == 0) "-"
                else format(signif(z$said / err, 3)), verdict))
  }
}
cat("\n-- tally --\n")
for (k in names(tally)) {
  t <- table(factor(tally[[k]], c("ok silent", "ok warn",
                                  "FALSE ALARM", "MISS")))
  cat(sprintf("  ss_extrapolate = %-5s  ok silent %d  ok warn %d",
              k, t[["ok silent"]], t[["ok warn"]]))
  cat(sprintf("  FALSE ALARM %d  MISS %d\n", t[["FALSE ALARM"]],
              t[["MISS"]]))
}
cat("\ndone\n")
