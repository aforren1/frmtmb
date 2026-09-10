# REVIEW of lane nss, attack 8: construct the case where the guarded
# thing is PRESENT and the guard is QUIET.
#
# rev-nss-06-backcompat.R found two of 26 schedules where the base
# commit warns and the lane's `ss_extrapolate = FALSE` arm does not.
# Both are structural, not numerical.
#
# The mechanism. On the FALSE arm `ode_run_in_finish()` reports
#
#     rel <- max(abs(ext[["y"]] - last)) / scale
#
# which is THE SIZE OF THE CORRECTION IT DECLINED TO MAKE. The three
# guards of `ode_ss_extrapolate()` all work by driving that correction
# to EXACTLY ZERO: `w` is 0 once `r` reaches 1, `rr` is 0 once `r` goes
# negative, and the damped denominator sends `r` to 0 when the
# differences are at the solver's noise. So on every state the guards
# were written for, the FALSE arm's `rel` is 0 and the check cannot
# fire, no matter how far the run-in is from its limit. The base commit
# compared the last two cycle-start states and had no such hole.
#
# The `stalled` flag does not rescue it: it is recorded INSIDE
#     if (is.finite(rel) && rel > ss_tol) { ... }
# so a `rel` of 0 suppresses it too.
#
# Second case: `n_ss` of 1. `keep` is filled from `j <- k - (n_ss - 3)`,
# so at n_ss = 1 only slot 3 is filled, `have[2] && have[3]` is FALSE,
# `ext` is NULL and `ode_run_in_finish()` returns before any check at
# all. The base commit compared `prev` with `y` at every n_ss >= 1.
#
# Script path: dev/rev-nss/rev-nss-07-warnopen.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

one_auc <- function(t, y, p) list(c(-p[2L] * y[1L],
                                    p[2L] * y[1L] - p[1L] * y[2L],
                                    y[2L]))
two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))

go <- function(tag, dyn, ns, pv, ii, out, n_ss, ext, tol = 1e-6) {
  args <- list(dyn, init = rep(list(0), ns),
               times = seq(0, ii, length.out = 9), parms = pv,
               events = data.frame(time = 0, state = 1L, value = 100,
                                   ii = ii, ss = TRUE),
               output = out, n_ss = n_ss, ss_tol = tol,
               atol = 1e-8, rtol = 1e-8)
  if (!is.null(ext)) args$ss_extrapolate <- ext
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(do.call(frm_ode, args),
             error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(e) { w <<- c(w, conditionMessage(e))
                            invokeRestart("muffleWarning") })
  hit <- any(grepl("run-in", w))
  cat(sprintf("  %-38s ext=%-5s n_ss=%3d  warns: %-5s\n", tag,
              if (is.null(ext)) "base" else ext, n_ss, hit))
  if (hit) cat("      ", substr(w[grepl("run-in", w)][1L], 1, 200),
               "\n")
  invisible(list(v = v, w = w, hit = hit))
}

cat("\n== A. an AUC compartment, which has NO steady state ==\n")
cat("The run-in's third state grows without bound every cycle. A",
    "warning\nis the only thing a user can act on.\n\n")
if (identical(arm, "ref")) {
  go("1 cmt oral + AUC, q12", one_auc, 3L, list(0.05, 1), 12, 3L, 20L,
     NULL)
  go("1 cmt oral + AUC, read the CENTRAL state", one_auc, 3L,
     list(0.05, 1), 12, 2L, 20L, NULL)
} else {
  go("1 cmt oral + AUC, q12", one_auc, 3L, list(0.05, 1), 12, 3L, 20L,
     FALSE)
  go("1 cmt oral + AUC, q12", one_auc, 3L, list(0.05, 1), 12, 3L, 20L,
     TRUE)
  go("1 cmt oral + AUC, read the CENTRAL state", one_auc, 3L,
     list(0.05, 1), 12, 2L, 20L, FALSE)
  go("1 cmt oral + AUC, read the CENTRAL state", one_auc, 3L,
     list(0.05, 1), 12, 2L, 20L, TRUE)
}

cat("\n== B. a short run-in, where `keep` is not full ==\n")
for (n in c(1L, 2L, 3L)) {
  if (identical(arm, "ref"))
    go("2 cmt oral 107 h, q24", two_oral, 3L,
       list(0.15, 0.3, 0.02, 1), 24, 2L, n, NULL)
  else for (e in c(FALSE, TRUE))
    go("2 cmt oral 107 h, q24", two_oral, 3L,
       list(0.15, 0.3, 0.02, 1), 24, 2L, n, e)
}

cat("\n== C. how far from the limit are those silent answers? ==\n")
lin <- function(ke, ka, ii, tt)
  as.numeric(frm_lincmt(parms = list(ke = ke, ka = ka, V = 10),
                        times = tt, ncmt = 1, depot = TRUE,
                        events = data.frame(time = 0, state = "depot",
                                            value = 100, ii = ii,
                                            addl = 0L, ss = TRUE)))
tt <- seq(0, 12, length.out = 9)
ref <- lin(0.05, 1, 12, tt)
for (e in if (identical(arm, "ref")) list(NULL) else list(FALSE, TRUE)) {
  z <- go("central state, 1 cmt oral + AUC", one_auc, 3L,
          list(0.05, 1), 12, 2L, 20L, e)
  if (!is.character(z$v))
    cat(sprintf("      central-state error against the exact limit:",
                ), sprintf(" %.3e\n",
                           max(abs(as.numeric(z$v) / 10 - ref)) /
                             max(abs(ref))))
}
cat("\ndone\n")
