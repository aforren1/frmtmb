# Lane `latent`, item 2.3, arm B: line the two models up exactly.
#
# Probe 2 found the 19.7 gap is not a likelihood difference: hmmTMB's
# `initial_state = "estimated"` estimates a SEPARATE initial
# distribution per sequence (20 of them here, 38 extra parameters),
# while frmtmb's `init = "estimated"` estimates ONE shared across
# sequences. Different models, so of course different likelihoods.
#
# This script pins the initial distribution to uniform on BOTH sides,
# which is the frmtmb spelling the Phase 0 scale row already uses and is
# also the truth the simulator draws from, and asks what is left.
#
#   Rscript dev/latent-2p3-hmmtmb-probe3.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(hmmTMB)
})
source("dev/latent-hmm-sim.R")

s <- hmm_sim(seed = 20260909L, ns = 20L, tl = 100L)
d <- s$d
K <- 3L
form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fam <- hmm(K = K, gaussian(), time = t, group = id, init = "uniform")
fit <- suppressWarnings(frm(form, family = fam, data = d))
ll_frm <- as.numeric(logLik(fit))
b <- unlist(fixef(fit))
sd_re <- sqrt(VarCorr(fit)[[1L]][1L, 1L])
nre <- nlevels(d$id)
cat("frm logLik:", format(ll_frm, digits = 12), "  sd(re):",
    format(sd_re, digits = 8), "  n re:", nre, "\n")

dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"
cat("\nformula matrix, with the reference column blanked:\n")
print(fmat)
hid <- MarkovChain$new(data = dh, n_states = K, formula = fmat,
                       initial_state = "estimated")
cat("ref:", hid$ref(), "\n")
cat("delta0:\n"); print(head(hid$delta0(), 3))
d0 <- hid$delta0()
cat("dim(delta0):", dim(d0), "\n")
cat("ref_delta0:", utils::head(hid$ref_delta0()), "\n")
cat("fixpar:\n"); print(hid$fixpar())

## Pin delta0 to uniform and FIX it, which is what frmtmb's
## init = "uniform" does.
d0[] <- 1 / K
hid$update_delta0(delta0 = d0)
fp <- hid$fixpar()
cat("\nfixpar after update_delta0:\n"); print(fp)
