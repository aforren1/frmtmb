# REVIEW of lane nss, attack 8, the sharp version.
#
# rev-nss-07-warnopen.R showed the shape. This one pins the size: on
# `ss_extrapolate = FALSE` the lane reports `rel` as the SIZE OF THE
# CORRECTION IT DID NOT MAKE, and every one of the three guards works by
# making that correction exactly zero. So on a state the guards fire on,
# the FALSE arm's check is a check that cannot fail, which is the exact
# defect the lane's own code comment says it fixed for the other arm.
#
# Three constructions, each with the value the user actually gets:
#   A. a state with no steady state, read AS THE OUTPUT;
#   B. a run-in that oscillates, so `r` is negative and `rr` is 0;
#   C. n_ss = 1, where `keep` is not full and no check runs at all.
#
# Script path: dev/rev-nss/rev-nss-08-failopen.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

run1 <- function(dyn, ns, pv, ii, out, n_ss, ext, tt) {
  args <- list(dyn, init = rep(list(0), ns), times = tt, parms = pv,
               events = data.frame(time = 0, state = 1L, value = 100,
                                   ii = ii, ss = TRUE),
               output = out, n_ss = n_ss, ss_tol = 1e-6,
               atol = 1e-10, rtol = 1e-10)
  if (!is.null(ext)) args$ss_extrapolate <- ext
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(do.call(frm_ode, args),
             error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(e) { w <<- c(w, conditionMessage(e))
                            invokeRestart("muffleWarning") })
  ww <- w[grepl("run-in", w)]
  list(v = if (is.character(v)) v else as.numeric(v),
       warned = length(ww) > 0L,
       rel = if (!length(ww)) NA_real_ else
         as.numeric(sub("^.*(?:moving by|is still) ([0-9.e+-]+) .*$",
                        "\\1", ww[[1L]])))
}

arms <- if (identical(arm, "ref")) list(base = NULL) else
  list(`ext=FALSE` = FALSE, `ext=TRUE` = TRUE)

report <- function(tag, dyn, ns, pv, ii, out, n_ss, tt, truth) {
  cat(sprintf("\n-- %s --\n", tag))
  for (an in names(arms)) {
    z <- run1(dyn, ns, pv, ii, out, n_ss, arms[[an]], tt)
    err <- if (is.character(z$v) || is.null(truth)) NA_real_ else
      max(abs(z$v - truth)) / max(abs(truth))
    cat(sprintf("   %-10s warned %-5s  reported rel %10s  ",
                an, z$warned,
                if (is.na(z$rel)) "-" else format(signif(z$rel, 3))))
    cat(sprintf("TRUE error of the returned value %s\n",
                if (is.na(err)) "-" else format(signif(err, 4))))
  }
}

## A. a state with no steady state, READ AS THE OUTPUT ----------------
# state 3 integrates the central concentration, so it grows by the same
# AUC every cycle for ever and has no limit at all.
one_auc <- function(t, y, p) list(c(-p[2L] * y[1L],
                                    p[2L] * y[1L] - p[1L] * y[2L],
                                    y[2L]))
cat("\n== A. an AUC state read as the output, q12, n_ss = 20 ==")
report("AUC state is the output", one_auc, 3L, list(0.05, 1), 12, 3L,
       20L, seq(0, 12, length.out = 9), NULL)

## B. a run-in that oscillates ----------------------------------------
# a damped oscillator dosed every ii: the cycle map has complex
# eigenvalues, the successive differences change sign, `r` goes
# negative and `rr` clamps it to 0, so the correction is exactly 0.
osc <- function(t, y, p) list(c(y[2L],
                                -p[1L] * p[1L] * y[1L] -
                                  2 * p[2L] * p[1L] * y[2L]))
cat("\n== B. a run-in that oscillates, so r < 0 and rr = 0 ==")
report("damped oscillator, ii = 1", osc, 2L, list(3.0, 0.05), 1, 1L,
       20L, seq(0, 1, length.out = 9), NULL)

## C. n_ss = 1, where `keep` is not full ------------------------------
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
lim <- as.numeric(frm_lincmt(
  parms = list(ke = 0.15, k12 = 0.3, k21 = 0.02, ka = 1, V = 10),
  times = seq(0, 24, length.out = 9), ncmt = 2, depot = TRUE,
  events = data.frame(time = 0, state = "depot", value = 100, ii = 24,
                      addl = 0L, ss = TRUE))) * 10
cat("\n== C. n_ss = 1 and 2, 2 cmt oral 107 h q24 ==")
for (n in c(1L, 2L))
  report(paste("n_ss =", n), two_oral, 3L, list(0.15, 0.3, 0.02, 1),
         24, 2L, n, seq(0, 24, length.out = 9), lim)

cat("\ndone\n")
