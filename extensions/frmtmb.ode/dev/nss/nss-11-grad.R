# Does the tail correction differentiate?
#
# Two questions, and they are different. (1) Does the tape compute the
# derivative of what the numeric path computes: AD against a central
# difference of frm_ode() itself. (2) Is that derivative the one the
# model wants: AD against a central difference of the EXACT steady
# state, which frm_lincmt() at n_ss = Inf supplies in closed form. The
# truncated run-in answers (1) and fails (2), which is what makes the
# defect a wrong answer during a fit rather than a slow one.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-11-grad.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
nss_report_env()

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}
II <- 24
TT <- c(0, 2, 6, 12, 23.9)
EV <- data.frame(time = 0, state = 1L, value = 100, ii = II, ss = TRUE)
EVL <- data.frame(time = 0, state = "depot", value = 100, ii = II,
                  addl = 0L, ss = TRUE)

f_ode <- function(th, ext, n = 20L) {
  sum(frm_ode(two_oral, init = list(0, 0, 0), times = TT,
              parms = list(exp(th[1L]), exp(th[2L]), exp(th[3L]),
                           exp(th[4L])),
              events = EV, output = 2L, n_ss = n, ss_tol = 1,
              ss_extrapolate = ext, atol = 1e-10, rtol = 1e-10))
}
f_exact <- function(th) {
  10 * sum(frm_lincmt(parms = list(ke = exp(th[1L]), k12 = exp(th[2L]),
                                   k21 = exp(th[3L]), ka = exp(th[4L]),
                                   V = 10),
                      times = TT, ncmt = 2, depot = TRUE, events = EVL))
}
fd <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    a <- x; b <- x
    a[j] <- a[j] + h; b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h)
  }, 0)
}

x0 <- log(c(0.15, 0.3, 0.02, 1.0))
cat("\ndesign: 2 cmt oral, terminal half-life 107 h, ii = 24,",
    "n_ss = 20,\natol = rtol = 1e-10, theta = log(ke, k12, k21, ka)\n\n")

for (ext in c(TRUE, FALSE)) {
  tp <- MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    f_ode(th, ext)
  }, x0)
  g_ad <- as.numeric(tp$jacfun()(x0))
  g_fd <- fd(function(th) f_ode(th, ext), x0)
  g_ex <- fd(f_exact, x0)
  v_ad <- tp(x0)
  v_ex <- f_exact(x0)
  cat(sprintf("ss_extrapolate = %s\n", ext))
  cat(sprintf("  value: tape %14.8f  exact %14.8f  rel %9.2e\n",
              v_ad, v_ex, abs(v_ad - v_ex) / abs(v_ex)))
  cat(sprintf("  AD vs its own numeric path, max rel  %9.2e\n",
              max(abs(g_ad - g_fd)) / max(abs(g_fd))))
  cat(sprintf("  AD vs the EXACT steady state, max rel %9.2e\n",
              max(abs(g_ad - g_ex)) / max(abs(g_ex))))
  cat("    d/dlog(k21): tape", format(g_ad[3L], digits = 9),
      " exact", format(g_ex[3L], digits = 9), "\n")
}

cat("\nthe same at n_ss = 5, where truncation is far worse:\n")
for (ext in c(TRUE, FALSE)) {
  tp <- MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    f_ode(th, ext, 5L)
  }, x0)
  g_ad <- as.numeric(tp$jacfun()(x0))
  g_fd <- fd(function(th) f_ode(th, ext, 5L), x0)
  g_ex <- fd(f_exact, x0)
  cat(sprintf("  ss_extrapolate = %-5s  AD vs own path %9.2e   ",
              ext, max(abs(g_ad - g_fd)) / max(abs(g_fd))))
  cat(sprintf("AD vs exact %9.2e\n",
              max(abs(g_ad - g_ex)) / max(abs(g_ex))))
}
