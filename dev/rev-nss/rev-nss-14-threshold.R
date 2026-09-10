# REVIEW of lane nss round 1: two sharp questions the broad sweep
# raised but could not answer.
#
# A. IS THERE A MISS? `dev/rev-nss/rev-nss-13.log` found none in 64
#    runs, but it also found the report understating the true error by
#    up to 96x on the default arm. A miss needs the reported number
#    below `ss_tol` while the true error is above it, so it lives where
#    the true error is between `ss_tol` and about 100 * `ss_tol`. None
#    of the 32 cases landed there by accident. This walks `n_ss` until
#    one does, on the two families with the worst understatement.
#
# B. DID THE LANE'S LOW-GATE PROBE REACH THE LOW GATE?
#    `dev/nss/nss-22-kink.R` sweeps `k21` so that
#    `exp(-lambda_z * ii)` = `ode_ss_rlow` = 0.05 and reports a gap
#    falling linearly, which it reads as C1 there. But at `n_ss` = 20 a
#    ratio of 0.05 leaves the run-in 0.05^20 = 1e-26 from its limit, so
#    the differences the run-in measures are far below the damping
#    floor and the ratio the function actually SEES is 0, not 0.05. The
#    gate is shut on both sides of that sweep. This prints the ratio
#    the function sees, and then repeats the sweep at an `n_ss` short
#    enough for the differences to survive.
#
# Script path: dev/rev-nss/rev-nss-14-threshold.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
rev_env()
SS_TOL <- 1e-6
RLOW <- frmtmb.ode:::ode_ss_rlow

one_oral <- function(t, y, p) list(c(-p[2L] * y[1L],
                                     p[2L] * y[1L] - p[1L] * y[2L]))
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))

go <- function(ncmt, pv, ii, n_ss, ext, tol = 1e-10) {
  dyn <- if (ncmt == 1L) one_oral else two_oral
  ns <- if (ncmt == 1L) 2L else 3L
  tt <- seq(0, ii, length.out = 25)
  w <- character(0)
  v <- withCallingHandlers(
    frm_ode(dyn, init = rep(list(0), ns), times = tt, parms = pv,
            events = data.frame(time = 0, state = 1L, value = 100,
                                ii = ii, ss = TRUE),
            output = 2L, n_ss = n_ss, ss_tol = SS_TOL,
            ss_extrapolate = ext, atol = tol, rtol = tol),
    warning = function(e) { w <<- c(w, conditionMessage(e))
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

## ---- A. walk n_ss into the miss window -----------------------------
fam <- list(
  list(tag = "big peripheral 100x, t1/2z 340 h", ncmt = 2L,
       pv = list(0.05, 0.50, 0.005, 1.0),
       P = list(ke = 0.05, k12 = 0.5, k21 = 0.005, ka = 1, V = 10),
       ii = 24, ns = c(20, 40, 60, 80, 100, 120, 150, 200, 260, 320)),
  list(tag = "ramp band, 1 cmt t1/2 1653 h", ncmt = 1L,
       pv = list(log(2) / 1653, 1.0),
       P = list(ke = log(2) / 1653, ka = 1, V = 10),
       ii = 24, ns = c(20, 60, 120, 200, 300, 400, 500, 650, 800,
                       1000)),
  list(tag = "big peripheral 400x, t1/2z 900 h", ncmt = 2L,
       pv = list(0.02, 0.80, 0.002, 1.0),
       P = list(ke = 0.02, k12 = 0.8, k21 = 0.002, ka = 1, V = 10),
       ii = 24, ns = c(60, 120, 200, 300, 400, 550, 700, 900)))
cat("\n== A. walking n_ss through the miss window,",
    "ss_extrapolate = TRUE ==\n")
cat("   a MISS is true error > ss_tol with the warning silent\n\n")
cat(sprintf("%-34s %5s %11s %11s %9s %8s %s\n", "case", "n_ss",
            "true error", "it said", "said/true", "warned", "verdict"))
nmiss <- 0L
for (f in fam) {
  ref <- lin(f$ncmt, f$P, f$ii)
  for (n in f$ns) {
    z <- go(f$ncmt, f$pv, f$ii, as.integer(n), TRUE)
    err <- max(abs(z$v - ref)) / max(abs(ref))
    bad <- err > SS_TOL
    verdict <- if (z$warned && bad) "ok warn"
               else if (!z$warned && !bad) "ok silent"
               else if (z$warned) "FALSE ALARM" else "MISS"
    if (verdict == "MISS") nmiss <- nmiss + 1L
    cat(sprintf("%-34s %5d %11.3e %11s %9s %8s %s\n", f$tag, n, err,
                if (is.na(z$said)) "-" else format(signif(z$said, 3)),
                if (is.na(z$said)) "-" else
                  format(signif(z$said / err, 3)), z$warned, verdict))
  }
  cat("\n")
}
cat("misses found:", nmiss, "\n")

## ---- B. does the lane's low-gate probe reach the low gate? ---------
cat("\n== B. the ratio ode_ss_extrapolate() actually sees ==\n")
# Rebuild the run-in exactly as ode_run_in() does, on the exact cycle
# map, and read the ratio the shipped function returns for it.
ext <- frmtmb.ode:::ode_ss_extrapolate
seen <- function(ke, ka, ii, n_ss, tol) {
  M <- matrix(c(-ka, 0, ka, -ke), 2, 2)
  A <- as.matrix(Matrix::expm(M * ii))
  e <- c(100, 0)
  y <- c(0, 0)
  keep <- vector("list", n_ss + 1L)
  keep[[1L]] <- y
  for (k in seq_len(n_ss)) {
    y <- as.numeric(A %*% (y + e)); keep[[k + 1L]] <- y
  }
  z <- ext(keep[[n_ss - 1L]], keep[[n_ss]], keep[[n_ss + 1L]], tol, tol)
  c(as.numeric(z$r)[2L],
    abs(keep[[n_ss]][2L] - keep[[n_ss - 1L]][2L]) / abs(y[2L]))
}
ke <- -log(RLOW) / 24
cat(sprintf("  1 cmt, ke %.5f, ka 1, ii 24: exp(-ke * ii) = %.5f\n",
            ke, exp(-ke * 24)))
cat(sprintf("%8s %16s %18s\n", "n_ss", "r the code sees",
            "|d1| / |y| central"))
for (n in c(3L, 4L, 5L, 8L, 12L, 20L)) {
  s <- seen(ke, 1, 24, n, 1e-10)
  cat(sprintf("%8d %16.6e %18.3e\n", n, s[[1L]], s[[2L]]))
}

cat("\n== B2. the same one-sided sweep at an n_ss the gate survives ==\n")
lamz <- function(ke) ke
k_at <- function(target, ii) -log(target) / ii
grad_at <- function(kev, n_ss) {
  tp <- MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    sum(frm_ode(one_oral, init = list(0, 0),
                times = c(0, 6, 12, 18, 23.9),
                parms = list(exp(th[1L]), 1.0),
                events = data.frame(time = 0, state = 1L, value = 100,
                                    ii = 24, ss = TRUE),
                output = 2L, n_ss = n_ss, ss_tol = Inf,
                ss_extrapolate = TRUE, atol = 1e-12, rtol = 1e-12))
  }, log(kev))
  as.numeric(tp$jacfun()(log(kev)))
}
k0 <- k_at(RLOW, 24)
for (n in c(3L, 20L)) {
  cat(sprintf("\n-- n_ss = %d, ke at r = %.3f --\n", n, RLOW))
  cat(sprintf("%10s %16s %16s %12s %10s\n", "eps", "d/dlke below",
              "d/dlke above", "gap", "gap/eps"))
  for (eps in 10^-(3:6)) {
    # r = exp(-ke * ii) falls as ke rises, so "below r" is ke ABOVE k0
    a <- grad_at(k0 * (1 + eps), n)
    b <- grad_at(k0 * (1 - eps), n)
    cat(sprintf("%10.0e %16.8f %16.8f %12.3e %10.2e\n", eps, a, b,
                abs(a - b), abs(a - b) / eps))
  }
}
cat("\ndone\n")
