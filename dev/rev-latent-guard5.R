# Reviewer: two loose ends.
#
#  Ia (redone). guard4's contrast was wrong: evaluating fit$obj$fn at a
#      WORSE point does not move `last.par.best`, because RTMB only
#      writes it when the value improves. The hazard the hmm-starts.R
#      header names is a refit that finds a BETTER point, so that is
#      what is constructed here, on the d4 fit whose better optimum is
#      8.099 units away.
#
#  IV. the cost instrument. The 175.7 s figure cannot be re-taken here:
#      it is 22 minutes on a quiet machine and two other lanes are
#      running. What CAN be checked is whether the instrument's own
#      control does its job, so the same script is run at a size where
#      three blocks cost a minute, on a machine that is NOT quiet. A
#      control that stays near 1.0 on a busy machine is not measuring
#      the machine; a control that moves is the instrument working.
#
#   Rscript dev/rev-latent-guard5.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

hr("Ia (redone). last.par.best, moved by a BETTER evaluation")
d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-2p3-repro81.rds")
form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit <- frm(form, family = fam, data = d$dat)
cat("cold logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
ci0 <- suppressWarnings(stats::confint(fit))
lp0 <- fit[["obj"]][["env"]][["last.par.best"]]

# the better optimum's outer parameters, reached by hmm_starts()
ms <- suppressWarnings(hmm_starts(fit, n = 6, jitter = 2, seed = 4101))
cat("best logLik:", format(as.numeric(logLik(ms$best)), digits = 12), "\n")
ci1 <- suppressWarnings(stats::confint(fit))
lp1 <- fit[["obj"]][["env"]][["last.par.best"]]
cat("AFTER hmm_starts(): the original's confint identical:",
    identical(ci0, ci1), "  last.par.best identical:",
    identical(lp0, lp1), "\n")

# and now what a REUSED objective does: evaluate the original's own
# tape at the better optimum's parameters, which IS an improvement
pbest <- as.numeric(ms$best[["opt"]][["par"]])
v <- fit[["obj"]][["fn"]](pbest)
lp2 <- fit[["obj"]][["env"]][["last.par.best"]]
ci2 <- suppressWarnings(stats::confint(fit))
cat("one fn() call at the BETTER point, value", format(-v, digits = 12),
    "\n")
cat("  last.par.best moved                :", !identical(lp0, lp2), "\n")
cat("  max |change| in last.par.best      :",
    format(max(abs(as.numeric(lp2) - as.numeric(lp0))), digits = 4), "\n")
cat("  the original fit's confint() moved :", !identical(ci0, ci2), "\n")
if (!identical(ci0, ci2)) {
  cat("  max |change| in the est column     :",
      format(max(abs(ci2[, "est"] - ci0[, "est"])), digits = 4), "\n")
  cat("  max |change| in the interval width :",
      format(max(abs((ci2[, "upr"] - ci2[, "lwr"]) -
                       (ci0[, "upr"] - ci0[, "lwr"]))), digits = 4), "\n")
}
cat("  logLik(fit) now reads              :",
    format(as.numeric(logLik(fit)), digits = 12), "\n")

hr("IV. the cost instrument, run at a size a review can afford")
source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-hmm-sim.R")
elapsed <- function(expr) {
  gc(FALSE); t0 <- Sys.time(); force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}
s <- hmm_sim(seed = 20260910L, ns = 12L, tl = 150L)
dd <- s$d
f2 <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fm <- hmm(K = 3, gaussian(), time = t, group = id, init = "uniform")
cat("rows:", nrow(dd), " sequences:", nlevels(dd$id), "\n")
t_fit <- elapsed(fitc <- suppressWarnings(frm(f2, family = fm, data = dd)))
t_sdr <- elapsed(suppressWarnings(stats::confint(fitc)))
cat("one frm() fit:", format(t_fit, digits = 4), "s, ",
    fitc$opt$evaluations[2L], "gradient calls\n")
cat("the deferred sdreport:", format(t_sdr, digits = 4), "s\n")
tt <- numeric(3); gr <- numeric(3)
for (b in 1:3) {
  ms2 <- NULL
  tt[b] <- elapsed(ms2 <- suppressWarnings(
    hmm_starts(fitc, n = 2, jitter = 2, seed = 4201)))
  gr[b] <- sum(ms2$table$gr_evals[-1L], na.rm = TRUE)
  cat(sprintf("block %d: %.2f s, %g gradient calls, best %.6f\n",
              b, tt[b], gr[b], as.numeric(logLik(ms2$best))))
}
cat("\nCONTROL, max / min       :", format(max(tt) / min(tt), digits = 4),
    "  (the lane reports 1.023 at 25 000 rows on a quiet machine)\n")
cat("gradient calls per block :", gr, "\n")
cat("identical across blocks  :", length(unique(gr)) == 1L,
    "  <- the load-independent check, which is the one to prefer\n")
cat("seconds per refit, min   :", format(min(tt) / 2, digits = 5), "\n")
cat("seconds per gradient call, refit:",
    format(min(tt) / gr[1L], digits = 4), "\n")
cat("seconds per gradient call, plain fit:",
    format(t_fit / fitc$opt$evaluations[2L], digits = 4), "\n")
cat("ratio                    :",
    format((min(tt) / gr[1L]) / (t_fit / fitc$opt$evaluations[2L]),
           digits = 4),
    "  (the lane's own numbers give 6.63 / 11.47 = 0.578)\n")
cat("other R processes on the machine right now:",
    length(system("tasklist /FI \"IMAGENAME eq Rscript.exe\" /NH",
                  intern = TRUE)), "lines of tasklist\n")
