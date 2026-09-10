# Lane `latent`, item 2.3: what scale should a jittered start use?
#
# `hmm_starts()` has to perturb the outer parameter vector by an amount
# a user can reason about. Two candidates, and this script prints what
# each one is worth on probe D4's cold start: the fit's own standard
# errors (a run-measured unit, but it can be tiny at a sharp local
# optimum) and the parameter magnitudes.
#
#   Rscript dev/latent-2p3-scale-probe.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})

d <- readRDS("dev/latent-2p3-repro81.rds")
dat <- d$dat
form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit <- frm(form, family = fam, data = dat)

cat("logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n\n")
cat("opt$par:\n"); print(fit$opt$par)
cat("\nnames(fit$estimates):", names(fit$estimates), "\n")
cat("\nvcov rows:\n"); print(rownames(vcov(fit)))
ci <- suppressWarnings(stats::confint(fit))
cat("\nconfint rows:\n"); print(rownames(ci))
cat("\nconfint:\n"); print(round(ci, 4))
se <- (ci[, 2L] - ci[, 1L]) / (2 * stats::qnorm(0.975))
cat("\nimplied SEs:\n"); print(round(se, 4))

# what the two modes differ by, in each unit
cold <- fit$opt$par
cat("\ncold-start optimum, outer vector:\n"); print(round(cold, 4))
