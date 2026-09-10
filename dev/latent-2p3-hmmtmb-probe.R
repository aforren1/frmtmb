# Lane `latent`, item 2.3, arm B: get the hmmTMB mapping right on a
# small design before paying for it at 50 x 500.
#
# Two conventions have to line up before the comparison means anything:
#
#  - the REFERENCE CELL. frmtmb references state 1 in every row of the
#    transition matrix; hmmTMB references the DIAGONAL by default, but
#    takes `ref =`, so `ref = rep(1, K)` makes the two parameterizations
#    the same one and `tr12` is then literally the same coefficient.
#  - the INITIAL DISTRIBUTION. frmtmb refuses `init = "stationary"` when
#    a transition carries a predictor, so both sides use the estimated
#    initial distribution here.
#
#   Rscript dev/latent-2p3-hmmtmb-probe.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(hmmTMB)
})
source("dev/latent-hmm-sim.R")

s <- hmm_sim(seed = 20260909L, ns = 20L, tl = 100L)
d <- s$d
cat("rows:", nrow(d), "\n")

form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fam <- hmm(K = 3, gaussian(), time = t, group = id, init = "estimated")
t0 <- Sys.time()
fit <- frm(form, family = fam, data = d)
cat("frm seconds:", as.numeric(difftime(Sys.time(), t0, units = "secs")),
    "\n")
cat("frm logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
print(unlist(fixef(fit)))
print(VarCorr(fit))

## hmmTMB. `state` must NOT be a column: hmmTMB reads it as known states
## (probe D3 of dev/hmm-feasibility.md).
dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
K <- 3L
# hmmTMB insists the DIAGONAL of the formula matrix is ".", which is
# also its reference cell, so `ref = rep(1, K)` is not available once a
# formula matrix is given. It does not matter: row 1's reference is
# state 1 under both conventions, so hmmTMB's 1 -> 2 coefficient IS
# frmtmb's `tr12`. Rows 2 and 3 are parameterized differently, and both
# are saturated intercept-only logits there, so the two span the same
# model and the comparison is made on the TPM rather than on those
# coefficients.
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"
cat("\nformula matrix hmmTMB is given:\n")
print(fmat)
hid <- try(MarkovChain$new(data = dh, n_states = K,
                           formula = fmat,
                           initial_state = "estimated"), silent = TRUE)
if (inherits(hid, "try-error")) {
  cat("MarkovChain$new failed:", attr(hid, "condition")$message, "\n")
} else {
  print(hid$formulas())
  cat("ref:", hid$ref(), "\n")
  obs <- Observation$new(data = dh, n_states = K, dists = list(y = "norm"),
                         par = list(y = list(
                           mean = sort(unname(unlist(fixef(fit))[
                             c("mu1.(Intercept)", "mu2.(Intercept)",
                               "mu3.(Intercept)")])),
                           sd = rep(exp(unname(unlist(fixef(fit))[
                             "sigma1.(Intercept)"])), K))))
  hm <- HMM$new(obs = obs, hid = hid)
  t0 <- Sys.time()
  hm$fit(silent = TRUE)
  cat("hmmTMB seconds:",
      as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")
  cat("hmmTMB llk:", format(hm$llk(), digits = 12), "\n")
  print(hm$par()$obspar[, , 1])
  cat("coeff_fe (hid):\n"); print(hid$coeff_fe())
  cat("sd_re (hid):\n"); print(hid$sd_re())
}
