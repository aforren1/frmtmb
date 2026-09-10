# Lane `latent`, item 2.3, arm B: explain the transition-matrix column.
#
# Arm B's `tpm_max_gap` is 0.039 to 0.180, where its log-likelihood
# identity is 1e-11 and its estimate identity 1e-04. That is not a
# disagreement, it is the harness comparing two different quantities,
# and this script says which:
#
#   frmtmb's number is the POPULATION transition matrix, built from the
#   tr{i}{j} fixed-effect intercepts with the random effect at zero;
#   hmmTMB's `par()$tpm[, , 1]` is the transition matrix at DATA ROW 1,
#   which belongs to sequence 1 and therefore carries sequence 1's own
#   random intercept.
#
# The check: rebuild frmtmb's row-1 logits with sequence 1's
# conditional mode added, and compare THAT with hmmTMB's `tpm[, , 1]`.
# The fixed-effect coefficients are compared directly at the same time,
# since row 1 references state 1 under both conventions.
#
#   Rscript dev/latent-2p3-tpmgap.R [seed]

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(hmmTMB)
})
source("dev/latent-hmm-sim.R")

args <- commandArgs(trailingOnly = TRUE)
SEED <- if (length(args)) as.integer(args[[1L]]) else 20260910L
K <- 3L
s <- hmm_sim(seed = SEED, ns = 50L, tl = 500L)
d <- s$d
fit <- suppressWarnings(suppressMessages(
  frm(bf(y ~ 1, tr12 ~ 1 + (1 | id)),
      family = hmm(K = K, gaussian(), time = t, group = id,
                   init = "uniform"), data = d)))
b <- unlist(fixef(fit))
cat("frm logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")

dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"
h0 <- suppressMessages(MarkovChain$new(data = dh, n_states = K,
                                       formula = fmat,
                                       initial_state = "estimated"))
ld <- h0$delta0(log = TRUE, as_matrix = FALSE)
hid <- suppressMessages(MarkovChain$new(
  data = dh, n_states = K, formula = fmat, initial_state = "estimated",
  fixpar = list(delta0 = stats::setNames(rep(NA_real_, length(ld)),
                                         rownames(ld)))))
obs <- Observation$new(data = dh, n_states = K, dists = list(y = "norm"),
                       par = list(y = list(mean = hmm_truth$mu,
                                           sd = rep(hmm_truth$sigma, K))))
hm <- HMM$new(obs = obs, hid = hid)
suppressWarnings(hm$fit(silent = TRUE))
cat("hmmTMB llk :", format(hm$llk(), digits = 12), "\n")
cat("relative gap:",
    format(abs(as.numeric(logLik(fit)) - hm$llk()) /
             abs(as.numeric(logLik(fit))), digits = 4), "\n\n")

## the fixed-effect coefficients of ROW 1, which reference state 1
## under both conventions
cf <- hid$coeff_fe()
cat("tr12  frm:", format(unname(b["tr12.(Intercept)"]), digits = 9),
    "  hmmTMB:", format(cf["S1>S2.(Intercept)", 1L], digits = 9), "\n")
cat("tr13  frm:", format(unname(b["tr13.(Intercept)"]), digits = 9),
    "  hmmTMB:", format(cf["S1>S3.(Intercept)", 1L], digits = 9), "\n")
se_tr <- sqrt(diag(vcov(fit)))[c("tr12_(Intercept)", "tr13_(Intercept)")]
cat("  the fit's own SEs on those two:",
    format(unname(se_tr), digits = 4), "\n")
cat("  gap as a fraction of one SE   :",
    format(max(abs(c(unname(b["tr12.(Intercept)"]) -
                       cf["S1>S2.(Intercept)", 1L],
                     unname(b["tr13.(Intercept)"]) -
                       cf["S1>S3.(Intercept)", 1L]))) / min(se_tr),
           digits = 4), "\n\n")

## the transition matrices, the two ways
eta_pop <- matrix(0, K, K - 1L)
for (i in seq_len(K)) {
  for (j in 2:K) {
    eta_pop[i, j - 1L] <- unname(b[paste0("tr", i, j, ".(Intercept)")])
  }
}
G_pop <- hmm_tpm(eta_pop)
u1 <- ranef(fit)[[1L]][1L, 1L]
eta_1 <- eta_pop
eta_1[1L, 1L] <- eta_1[1L, 1L] + u1
G_1 <- hmm_tpm(eta_1)
G_h <- matrix(as.numeric(hm$par()$tpm[, , 1]), K, K)

cat("sequence 1's conditional mode on tr12:", format(u1, digits = 6),
    "\n\n")
cat("frmtmb POPULATION tpm:\n"); print(round(G_pop, 6))
cat("\nfrmtmb tpm WITH sequence 1's mode on tr12:\n")
print(round(G_1, 6))
cat("\nhmmTMB par()$tpm[, , 1]:\n"); print(round(G_h, 6))
cat("\nmax |population - hmmTMB[, , 1]|      :",
    format(max(abs(G_pop - G_h)), digits = 4), "\n")
cat("max |with-sequence-1 - hmmTMB[, , 1]| :",
    format(max(abs(G_1 - G_h)), digits = 4), "\n")
