# Lane `latent`, item 2.3: does frm_allfit() find probe D4's better
# optimum?
#
# `?hmm` used to send a reader to frm_allfit() for multimodality.
# frm_allfit() varies the OPTIMIZER from one fixed start, so it should
# not, but "should not" is not a measurement. This runs it on the same
# fit hmm_starts() recovers from.
#
#   Rscript dev/latent-2p3-allfit.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})

d <- readRDS("dev/latent-2p3-repro81.rds")
form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit <- frm(form, family = fam, data = d$dat)
cat("cold start logLik:", format(as.numeric(logLik(fit)), digits = 12),
    "\n")
cat("known optimum    : -1087.99646521\n\n")
af <- frm_allfit(fit)
print(af)
