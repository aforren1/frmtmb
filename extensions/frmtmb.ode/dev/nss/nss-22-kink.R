# NIT 4: is the default objective C1 now?
#
# The reviewer's design, from `dev/rev-nss/rev-nss-09-kink.R`,
# reproduced here so it runs against the lane library. Two-compartment
# oral, ke 0.15, k12 0.3, ka 1.0, ii 24, n_ss 20, atol = rtol = 1e-12,
# sweeping k21 across the value that puts exp(-lambda_z * ii) exactly
# at a target r. The test is not "is there a jump in the objective",
# which there never was, but whether the two ONE-SIDED derivative
# limits converge on each other as the bracket tightens.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-22-kink.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
nss_report_env()

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
k21_at <- function(rr) {
  L <- -log(rr) / II
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

sweep_at <- function(target, ext) {
  k0 <- k21_at(target)
  cat(sprintf("\n-- r = %g, k21 = %s, terminal half-life %s h,",
              target, format(k0, digits = 7),
              format(round(log(2) / lamz(k0), 1))))
  cat(sprintf(" ss_extrapolate = %s --\n", ext))
  cat(sprintf("%10s %16s %16s %14s\n", "eps", "d/dlk21 below",
              "d/dlk21 above", "gap"))
  for (eps in c(1e-3, 1e-4, 1e-5, 1e-6)) {
    lo <- grad(log(k0 * (1 - eps)), ext)
    hi <- grad(log(k0 * (1 + eps)), ext)
    cat(sprintf("%10.0e %16.6f %16.6f %14.4g\n", eps, lo, hi,
                abs(hi - lo)))
  }
}
cat("\n== across the stand-down point, r = rcap ==\n")
sweep_at(0.9875, TRUE)
sweep_at(0.9875, FALSE)
cat("\n== across r = 1, where the correction reaches zero ==\n")
sweep_at(1 - 1e-9, TRUE)
cat("\n== across the low gate at r = ode_ss_rlow ==\n")
sweep_at(0.05, TRUE)

cat("\n== what the smoothing costs: the largest factor the correction",
    "can apply ==\n")
# the shipped function itself, through the construction of
# nss-26-corners.R, so this cannot drift from what is installed
fac <- function(r) vapply(r, function(x) {
  if (x == 0) return(0)
  z <- frmtmb.ode:::ode_ss_extrapolate(0, 1 / x, 1 / x + 1, 1e-300,
                                       1e-300)
  as.numeric(z[["y"]]) - (1 / x + 1)
}, 0)
old <- function(r, rcap = 0.99) {
  w <- pmin(pmax((1 - r) / (1 - rcap), 0), 1)
  rr <- pmin(pmax(r, 0), rcap)
  w * rr / (1 - rr)
}
rg <- seq(0, 1.2, length.out = 4001)
cat(sprintf("  shipped now: max %.3f at r = %.5f\n", max(fac(rg)),
            rg[which.max(fac(rg))]))
cat(sprintf("  punch-0 form: max %.3f at r = %.5f\n", max(old(rg)),
            rg[which.max(old(rg))]))
cat(sprintf("  the two agree to %.3g for r <= 0.98\n",
            max(abs(fac(rg[rg <= 0.98]) - old(rg[rg <= 0.98])))))
