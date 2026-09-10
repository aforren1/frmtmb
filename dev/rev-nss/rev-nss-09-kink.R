# REVIEW of lane nss, a property the lane did not measure: the shipped
# correction makes the objective NON-SMOOTH in the parameters.
#
# `ode_ss_extrapolate()` builds min() and max() out of abs():
#
#     lo(x, a) = (x + a + |x - a|) / 2      max
#     hi(x, a) = (x + a - |x - a|) / 2      min
#     w  = max(min((1 - r) / (1 - 0.99), 1), 0)
#     rr = max(min(r, 0.99), 0)
#
# `abs()` is taped, which is why this works at all, and its derivative
# jumps at 0. So the correction `d2 * w * rr / (1 - rr)` is continuous
# in `r` but its derivative is not, at r = 0.99 and at r = 0 and at
# r = 1. Left and right of r = 0.99 the derivative of w * rr / (1 - rr)
# is +1 / (1 - r)^2 = +1e4 and -100 * 0.99 / 0.01 = -9900, a jump of
# about 2e4 times the last difference.
#
# This matters because the lane rejected recomputing `n_ss` inside the
# optimizer on exactly this ground: "the objective would then be
# discontinuous in the parameters wherever the chosen n_ss changed,
# which the Laplace approximation's inner Newton solve cannot carry".
# The shipped design is C0 rather than discontinuous, which is better,
# but it is not C1, and quasi-Newton and the inner Newton solve both
# assume C2.
#
# Is r = 0.99 reachable? r = exp(-lambda_z * ii), so r = 0.99 at
# ii = 24 is a terminal half-life of 1655 h. A user does not usually
# write that down, but an OPTIMIZER walks there: in
# dev/rev-nss/rev-nss-05.log seed 103 the exact-limit fit ran k21 to
# 3e-14, which is lambda_z at 1e-13 and r at 1 to fourteen figures. The
# path from a start to that boundary crosses 0.99.
#
# Script path: dev/rev-nss/rev-nss-09-kink.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
rev_env()

two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
II <- 24
TT <- c(0, 2, 6, 12, 23.9)
EV <- data.frame(time = 0, state = 1L, value = 100, ii = II, ss = TRUE)
KE <- 0.15; K12 <- 0.3; KA <- 1.0

lamz <- function(k21) {
  b <- KE + K12 + k21
  (b - sqrt(b * b - 4 * KE * k21)) / 2
}
# k21 that puts exp(-lambda_z * ii) exactly at the target r
k21_at <- function(rr) {
  L <- -log(rr) / II
  # lambda_z solves L^2 - (ke + k12 + k21) L + ke k21 = 0
  (L * L - (KE + K12) * L) / (L - KE)
}

f_ode <- function(th, ext, n = 20L)
  sum(frm_ode(two_oral, init = list(0, 0, 0), times = TT,
              parms = list(KE, K12, exp(th[1L]), KA), events = EV,
              output = 2L, n_ss = n, ss_tol = Inf,
              ss_extrapolate = ext, atol = 1e-12, rtol = 1e-12))
grad <- function(th, ext) {
  tp <- MakeTape(function(x) {
    "c" <- RTMB::ADoverload("c")
    f_ode(x, ext)
  }, th)
  as.numeric(tp$jacfun()(th))
}

cat("\n== the objective and its AD gradient across r = 0.99 ==\n")
cat("k21 at r = 0.99 is", format(k21_at(0.99)), ", terminal half-life",
    format(round(log(2) / lamz(k21_at(0.99)), 1)), "h\n\n")
k0 <- k21_at(0.99)
cat(sprintf("%12s %10s %14s %14s %14s\n", "k21", "r", "objective",
            "AD d/dlk21", "fd d/dlk21"))
prev <- NULL
for (m in c(0.7, 0.85, 0.95, 0.99, 0.999, 1.001, 1.01, 1.05, 1.15,
            1.4)) {
  k <- k0 * m
  th <- log(k)
  r <- exp(-lamz(k) * II)
  v <- f_ode(th, TRUE)
  g <- grad(th, TRUE)
  h <- 1e-6
  gf <- (f_ode(th + h, TRUE) - f_ode(th - h, TRUE)) / (2 * h)
  cat(sprintf("%12.3e %10.6f %14.7f %14.5f %14.5f\n", k, r, v, g, gf))
}

cat("\n== the same sweep with ss_extrapolate = FALSE, which is smooth ==\n")
cat(sprintf("%12s %10s %14s %14s\n", "k21", "r", "objective",
            "AD d/dlk21"))
for (m in c(0.7, 0.95, 0.999, 1.001, 1.05, 1.4)) {
  k <- k0 * m
  th <- log(k)
  cat(sprintf("%12.3e %10.6f %14.7f %14.5f\n", k, exp(-lamz(k) * II),
              f_ode(th, FALSE), grad(th, FALSE)))
}

cat("\n== the size of the jump, bracketing r = 0.99 tightly ==\n")
for (eps in c(1e-3, 1e-4, 1e-5)) {
  gl <- grad(log(k21_at(0.99 - eps)), TRUE)
  gr <- grad(log(k21_at(0.99 + eps)), TRUE)
  vl <- f_ode(log(k21_at(0.99 - eps)), TRUE)
  vr <- f_ode(log(k21_at(0.99 + eps)), TRUE)
  cat(sprintf("  eps %6.0e  r=%.5f g=%12.4f | r=%.5f g=%12.4f",
              eps, 0.99 - eps, gl, 0.99 + eps, gr))
  cat(sprintf("   objective gap %.3e\n", abs(vr - vl)))
}
cat("\ndone\n")
