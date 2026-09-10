# Reviewer, item 2.3 claim 1, last corner: the design rationale, and the
# usability of what comes back.
#
#   I. hmm-starts.R's header says each refit re-evaluates the call so
#      that RTMB's `last.par.best` is not moved and the ORIGINAL fit's
#      standard errors survive. That is a claim about behaviour and the
#      suite does not test it. Tested here.
#  II. $best is a fit a user will call summary(), confint() and vcov()
#      on. It is produced with `se = FALSE`; check those work.
# III. the two tolerances inside the function, side by side: the modes
#      table merges optima at grad_tol^2 relative, and the sentence
#      printed above it fires at any positive difference at all.
#
#   Rscript dev/rev-latent-guard4.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}
dd <- small_data()
fit <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)

hr("I. does hmm_starts() move the original fit's standard errors?")
ci_before <- suppressWarnings(stats::confint(fit))
v_before <- suppressWarnings(stats::vcov(fit))
lp_before <- fit[["obj"]][["env"]][["last.par.best"]]
ms <- suppressWarnings(hmm_starts(fit, n = 4, jitter = 2, seed = 55L))
ci_after <- suppressWarnings(stats::confint(fit))
v_after <- suppressWarnings(stats::vcov(fit))
lp_after <- fit[["obj"]][["env"]][["last.par.best"]]
cat("confint() bitwise identical before and after:",
    identical(ci_before, ci_after), "\n")
cat("vcov()    bitwise identical before and after:",
    identical(v_before, v_after), "\n")
cat("obj$env$last.par.best identical             :",
    identical(lp_before, lp_after), "\n")
cat("logLik(fit) unchanged                       :",
    identical(as.numeric(stats::logLik(fit)), ms$original_logLik), "\n")

hr("Ia. and what a REUSED objective would have done, for contrast")
lp0 <- fit[["obj"]][["env"]][["last.par.best"]]
junk <- fit[["obj"]][["fn"]](as.numeric(fit[["opt"]][["par"]]) + 0.5)
lp1 <- fit[["obj"]][["env"]][["last.par.best"]]
cat("one fn() call away from the optimum moved last.par.best:",
    !identical(lp0, lp1), "\n")
cat("  max |change| in last.par.best:",
    format(max(abs(as.numeric(lp1) - as.numeric(lp0))), digits = 4), "\n")
cat("  so the header's reason for re-evaluating the call is a real\n")
cat("  hazard, not a hypothetical one.\n")

hr("II. is $best a usable fit?")
dd2 <- small_data(seed = 4507L)
f2 <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
          data = dd2)
ms2 <- suppressWarnings(hmm_starts(f2, n = 4, jitter = 2, seed = 56L))
cat("best is a refit rather than the incumbent:",
    !identical(as.numeric(logLik(ms2$best)),
               as.numeric(logLik(f2))), "\n")
for (nm in c("summary", "confint", "vcov", "fixef", "logLik")) {
  r <- tryCatch({
    suppressWarnings(do.call(nm, list(ms2$best)))
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-8s on $best : %s\n", nm,
              if (identical(r, "ok")) "ok" else r))
}
cat("  frmtmb::diagnose() on $best:",
    tryCatch({
      suppressWarnings(frmtmb::diagnose(ms2$best, quiet = TRUE)); "ok"
    }, error = function(e) paste("ERROR:", conditionMessage(e))), "\n")

hr("III. the two tolerances, side by side")
cat("grad_tol default          :", 1e-3, "\n")
cat("mode_tol = grad_tol^2     :", 1e-6, " relative\n")
llx <- as.numeric(logLik(fit))
cat("on a fit at logLik", format(llx, digits = 8),
    "two optima are MERGED when they\n  differ by less than",
    format(1e-6 * abs(llx), digits = 4), "log-likelihood units.\n")
cat("the sentence 'the original fit found a local optimum' fires when\n")
cat("  the best refit is ANY amount above the incumbent: the test is\n")
cat("  `gap > 0` in print.frmtmb_hmm_starts(), and the swap that makes\n")
cat("  the gap is `lli > logLik(best)` in hmm_starts(). Neither has a\n")
cat("  tolerance.\n")
cat("measured gaps that triggered it on unimodal fits (guard3, case F):\n")
cat("  6.0e-11 to 1.5e-09, against a merge threshold of",
    format(1e-6 * abs(llx), digits = 3), "\n")
cat("ratio, threshold over the largest false trigger:",
    format(1e-6 * abs(llx) / 1.4642e-09, digits = 5), "\n")
