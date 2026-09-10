# The bars the new tests assert, computed with the shipped guard and
# with each piece of it ABSENT, so that each test is known to
# discriminate rather than assumed to.
#
# The systems are the ones the tests use, and the run-in is replayed
# here rather than inside frm_ode() so that a variant rule can be
# substituted for the shipped one.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-14-absent.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
cat("frmtmb.ode from", dirname(find.package("frmtmb.ode")), "\n")

ATOL <- 1e-10
RTOL <- 1e-10
run_in <- function(dyn, p, ns, ii, amt, n) {
  y <- numeric(ns)
  ys <- matrix(0, n + 1L, ns)
  for (k in seq_len(n)) {
    y[1L] <- y[1L] + amt
    s <- deSolve::ode(y, c(0, ii), dyn, p, method = "lsoda",
                      atol = ATOL, rtol = RTOL)
    y <- as.numeric(s[2L, 1L + seq_len(ns)])
    ys[k + 1L, ] <- y
  }
  ys
}
lo <- function(x, a) (x + a + abs(x - a)) / 2
hi <- function(x, a) (x + a - abs(x - a)) / 2

# the shipped rule, then each piece of it removed in turn
variants <- list(
  truncated = function(ya, yb, yc) yc,
  shipped = function(ya, yb, yc) {
    d1 <- yb - ya; d2 <- yc - yb
    del <- (4 * (ATOL + RTOL * abs(yc)))^2
    r <- (d1 * d2) / (d1 * d1 + del)
    w <- lo(hi((1 - r) / (1 - 0.99), 1), 0)
    rr <- lo(hi(r, 0.99), 0)
    yc + d2 * w * rr / (1 - rr)
  },
  no_perstate = function(ya, yb, yc) {
    d1 <- yb - ya; d2 <- yc - yb
    del <- sum((4 * (ATOL + RTOL * abs(yc)))^2)
    r <- sum(d1 * d2) / (sum(d1 * d1) + del)
    w <- lo(hi((1 - r) / (1 - 0.99), 1), 0)
    rr <- lo(hi(r, 0.99), 0)
    yc + d2 * w * rr / (1 - rr)
  },
  no_damping = function(ya, yb, yc) {
    d1 <- yb - ya; d2 <- yc - yb
    r <- (d1 * d2) / (d1 * d1)
    w <- lo(hi((1 - r) / (1 - 0.99), 1), 0)
    rr <- lo(hi(r, 0.99), 0)
    yc + d2 * w * rr / (1 - rr)
  },
  no_ramp = function(ya, yb, yc) {
    d1 <- yb - ya; d2 <- yc - yb
    del <- (4 * (ATOL + RTOL * abs(yc)))^2
    r <- (d1 * d2) / (d1 * d1 + del)
    rr <- lo(hi(r, 0.99), 0)
    yc + d2 * rr / (1 - rr)
  })

report <- function(tag, dyn, p, ns, ii, keep, nref) {
  ys <- run_in(dyn, p, ns, ii, 100, nref)
  ref <- ys[nref + 1L, keep]
  cat(sprintf("\n%s (reading state %s, reference n_ss = %d)\n", tag,
              paste(keep, collapse = ","), nref))
  for (v in names(variants)) {
    y <- variants[[v]](ys[19L, ], ys[20L, ], ys[21L, ])
    e <- max(abs(y[keep] - ref)) / max(abs(ref))
    cat(sprintf("  %-12s error %11.3e\n", v, e))
  }
}

# T2: a slow one-compartment oral system, where the tail is the defect
slow1 <- function(t, y, p) list(c(-p$ka * y[1L],
                                  p$ka * y[1L] - p$ke * y[2L]))
report("T2 slow 1 cmt oral, ke*ii = 0.12", slow1,
       list(ka = 1, ke = 0.01), 2L, 12, 2L, 4000L)

# T4: the same system with an AUC state riding along, which has no
# steady state at all
slow1auc <- function(t, y, p) list(c(-p$ka * y[1L],
                                     p$ka * y[1L] - p$ke * y[2L],
                                     y[2L]))
report("T4 the same plus an AUC state", slow1auc,
       list(ka = 1, ke = 0.01), 3L, 12, 2L, 4000L)

# T5: a state the dosing never reaches, so its difference is exactly
# zero every cycle and an undamped ratio is 0 / 0
dead <- function(t, y, p) list(c(-p$ka * y[1L],
                                 p$ka * y[1L] - p$ke * y[2L],
                                 0 * y[3L]))
report("T5 an unreachable state, difference exactly 0", dead,
       list(ka = 1, ke = 0.01), 3L, 12, 2L, 4000L)

# T6: a component that is GROWING between cycles, which a cap reads as
# "nearly stopped" and the ramp reads as "do not extrapolate this one"
ATOL <- 1e-8
RTOL <- 1e-8
osc <- function(t, y, p) list(c(y[2L], -p$w2 * y[1L] - p$z * y[2L]))
report("T6 oscillator, ii = 1, atol = rtol = 1e-8", osc,
       list(w2 = 1, z = 0.1), 2L, 1, 1L, 8000L)
