# REVIEW of lane nss: reproduce section 4's table, and measure the one
# thing it does not say.
#
# The correction multiplies the run-in's LAST DIFFERENCE by r / (1 - r),
# which is 99 at the cap. That difference carries the integrator's own
# error, so the answer the correction returns carries up to 99 times
# it. The lane writes "the extrapolated column is the INTEGRATOR's
# tolerance, not the run-in's", and ?frm_ode repeats it. That is right
# where r is small and wrong by about 1 / (1 - r) where it is not, so
# the amplification is reported here beside the error.
#
# Both are measured against `frm_lincmt()` at its default n_ss = Inf,
# which is the exact limit, over 25 points of one dosing interval, at
# the shipped n_ss = 20 and three tolerance settings.
#
# Script path: dev/rev-nss/rev-nss-11-amplify.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}

cases <- list(
  list(t = "23.2 h", ii = 8,  p = c(0.2, 0.4, 0.1, 1.1)),
  list(t = "107.1 h", ii = 24, p = c(0.15, 0.3, 0.02, 1.0)),
  list(t = "107.1 h", ii = 12, p = c(0.15, 0.3, 0.02, 1.0)),
  list(t = "264.6 h", ii = 24, p = c(0.1, 0.2, 0.008, 1.0)),
  list(t = "669.9 h", ii = 24, p = c(0.08, 0.15, 0.003, 1.0)),
  list(t = "1655 h",  ii = 24, p = c(0.15, 0.3, 0.001258637, 1.0)),
  list(t = "3.5 h",   ii = 24, p = c(0.2, 0.4, 0.5, 1.0)))

cat(sprintf("\n%-9s %4s %9s %9s %11s %11s %9s %9s\n", "t1/2z", "ii",
            "lz*ii", "r", "truncated", "extrapolated", "gain",
            "err/tol"))
for (tol in c(1e-6, 1e-8, 1e-10)) {
  cat(sprintf("\n-- atol = rtol = %g --\n", tol))
  for (cs in cases) {
    lz <- lamz(cs$p[1L], cs$p[2L], cs$p[3L])
    tt <- seq(0, cs$ii, length.out = 25)
    ev <- data.frame(time = 0, state = 1L, value = 100, ii = cs$ii,
                     ss = TRUE)
    evl <- data.frame(time = 0, state = "depot", value = 100,
                      ii = cs$ii, addl = 0L, ss = TRUE)
    ref <- as.numeric(frm_lincmt(
      parms = list(ke = cs$p[1L], k12 = cs$p[2L], k21 = cs$p[3L],
                   ka = cs$p[4L], V = 10),
      times = tt, ncmt = 2, depot = TRUE, events = evl))
    got <- function(ext) suppressWarnings(as.numeric(
      frm_ode(two_oral, init = list(0, 0, 0), times = tt,
              parms = as.list(cs$p), events = ev, output = 2L,
              n_ss = 20L, ss_tol = Inf, ss_extrapolate = ext,
              atol = tol, rtol = tol))) / 10
    et <- max(abs(got(FALSE) - ref)) / max(abs(ref))
    ee <- max(abs(got(TRUE) - ref)) / max(abs(ref))
    cat(sprintf("%-9s %4d %9.4f %9.6f %11.3e %11.3e %9.1e %9.1f\n",
                cs$t, cs$ii, lz * cs$ii, exp(-lz * cs$ii), et, ee,
                et / ee, ee / tol))
  }
}
cat("\nerr/tol is the extrapolated error in units of the integrator's",
    "own\ntolerance. 1 supports \"the extrapolated column is the",
    "integrator's\ntolerance\"; 1 / (1 - r) is what the geometric sum",
    "predicts.\n")
cat("\ndone\n")
