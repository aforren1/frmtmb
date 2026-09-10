# Reviewer, item 2.3, attack 4: the two hmmTMB traps.
#
#   T1. `initial_state = "estimated"` fits one initial distribution PER
#       SEQUENCE where frmtmb fits one shared. Claim: 38 extra
#       parameters and 21.2 log-likelihood units on 20 sequences, with
#       nothing in either output saying so. Verified against the CORRECT
#       explanation, not the retraction: the retracted one, a missing
#       normalizer on hmmTMB's random-effect penalty, is also priced
#       here so the reader can see it does not close the gap.
#   T2. arm B's `tpm_max_gap` of 0.18 beside a 1e-11 identity is the
#       harness comparing frmtmb's POPULATION transition matrix with
#       hmmTMB's row-1 matrix, which carries sequence 1's own random
#       intercept. Checked at 20 x 100 rather than at the lane's
#       25 000 rows, because the claim is about which quantity is being
#       compared and that does not depend on the sample size.
#
#   Rscript dev/rev-latent-hmmtmb.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(hmmTMB)
})
source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-hmm-sim.R")
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

K <- 3L
NS <- 20L
s <- hmm_sim(seed = 20260909L, ns = NS, tl = 100L)
d <- s$d
dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"

mk_obs <- function() {
  Observation$new(data = dh, n_states = K, dists = list(y = "norm"),
                  par = list(y = list(mean = c(-2, 0, 3),
                                      sd = rep(0.7, K))))
}

hr("T1a. hmmTMB initial_state = 'estimated' against frmtmb's")
fit_e <- suppressWarnings(suppressMessages(
  frm(form, family = hmm(K = K, gaussian(), time = t, group = id,
                         init = "estimated"), data = d)))
ll_e <- as.numeric(logLik(fit_e))
hid_e <- suppressMessages(MarkovChain$new(data = dh, n_states = K,
                                          formula = fmat,
                                          initial_state = "estimated"))
hm_e <- HMM$new(obs = mk_obs(), hid = hid_e)
suppressWarnings(hm_e$fit(silent = TRUE))
ld <- hid_e$delta0(log = TRUE, as_matrix = FALSE)
cat("sequences                          :", NS, "\n")
cat("hmmTMB free delta0 entries         :", length(ld),
    "  ( (K-1) x n_seq =", (K - 1L) * NS, ")\n")
cat("frmtmb free initial-distribution df:", K - 1L, "\n")
cat("extra parameters hmmTMB carries    :", length(ld) - (K - 1L),
    "  claim 38\n")
cat("frmtmb logLik  :", format(ll_e, digits = 12), "\n")
cat("hmmTMB llk()   :", format(hm_e$llk(), digits = 12), "\n")
cat("hmmTMB minus frmtmb:", format(hm_e$llk() - ll_e, digits = 6),
    "  claim 21.2\n")
cat("does either output say the two models differ? checked by hand:\n")
cat("  frmtmb df attr(logLik)  :", attr(logLik(fit_e), "df"), "\n")
cat("  hmmTMB edf()            :", format(hm_e$edf(), digits = 8), "\n")

hr("T1b. the RETRACTED explanation, priced")
sd_h <- as.numeric(hid_e$sd_re()[1L, 1L])
nc <- -(NS / 2) * log(2 * pi * sd_h^2)
cat("hmmTMB sd_re                       :", format(sd_h, digits = 8), "\n")
cat("-(k/2) log(2 pi sd^2), k =", NS, "     :", format(nc, digits = 8),
    "  findings say 5.66 at the lane's own draw\n")
cat("gap                                :",
    format(hm_e$llk() - ll_e, digits = 8), "\n")
cat("residual after the normalizer      :",
    format(hm_e$llk() - ll_e - abs(nc), digits = 8),
    " <- not zero, so the normalizer is not the explanation\n")

hr("T1c. the CORRECT explanation: pin delta0 on both sides")
fit_u <- suppressWarnings(suppressMessages(
  frm(form, family = hmm(K = K, gaussian(), time = t, group = id,
                         init = "uniform"), data = d)))
ll_u <- as.numeric(logLik(fit_u))
fp <- stats::setNames(rep(NA_real_, length(ld)), rownames(ld))
hid_u <- suppressMessages(MarkovChain$new(data = dh, n_states = K,
                                          formula = fmat,
                                          initial_state = "estimated",
                                          fixpar = list(delta0 = fp)))
hm_u <- HMM$new(obs = mk_obs(), hid = hid_u)
suppressWarnings(hm_u$fit(silent = TRUE))
cat("frmtmb init = 'uniform'            :", format(ll_u, digits = 14), "\n")
cat("hmmTMB delta0 pinned at uniform    :", format(hm_u$llk(), digits = 14),
    "\n")
cat("absolute difference                :",
    format(abs(ll_u - hm_u$llk()), digits = 4), "\n")
cat("relative difference                :",
    format(abs(ll_u - hm_u$llk()) / abs(ll_u), digits = 4),
    "  claim 1e-11 or better\n")
cat("is delta0 actually uniform there?  :\n")
print(round(hid_u$delta0()[1:2, , drop = FALSE], 8))

hr("T2. the tpm column: population, or sequence 1's own matrix?")
b <- unlist(fixef(fit_u))
eta_pop <- matrix(0, K, K - 1L)
for (i in seq_len(K)) {
  for (j in 2:K) {
    eta_pop[i, j - 1L] <- unname(b[paste0("tr", i, j, ".(Intercept)")])
  }
}
G_pop <- hmm_tpm(eta_pop)
u1 <- ranef(fit_u)[[1L]][1L, 1L]
eta_1 <- eta_pop
eta_1[1L, 1L] <- eta_1[1L, 1L] + u1
G_1 <- hmm_tpm(eta_1)
G_h <- matrix(as.numeric(hm_u$par()$tpm[, , 1]), K, K)
cat("sequence 1's conditional mode on tr12:", format(u1, digits = 6), "\n")
cat("\nfrmtmb POPULATION tpm:\n"); print(round(G_pop, 6))
cat("frmtmb tpm + sequence 1's mode:\n"); print(round(G_1, 6))
cat("hmmTMB par()$tpm[, , 1]:\n"); print(round(G_h, 6))
cat("\nmax |population   - hmmTMB[, , 1]| :",
    format(max(abs(G_pop - G_h)), digits = 4), "\n")
cat("max |with mode u1 - hmmTMB[, , 1]| :",
    format(max(abs(G_1 - G_h)), digits = 4),
    " <- the harness artifact, confirmed\n")
cat("ratio of the two                   :",
    format(max(abs(G_pop - G_h)) / max(abs(G_1 - G_h)), digits = 4), "\n")
## and the check that it really is SEQUENCE 1 and not just any mode
u <- ranef(fit_u)[[1L]][, 1L]
cat("\nconditional modes, first five      :",
    paste(format(u[1:5], digits = 4), collapse = " "), "\n")
gaps <- vapply(seq_along(u), function(i) {
  e <- eta_pop; e[1L, 1L] <- e[1L, 1L] + u[i]
  max(abs(hmm_tpm(e) - G_h))
}, 1)
cat("which sequence's mode closes it    :", which.min(gaps),
    " (gap", format(min(gaps), digits = 3), ")\n")
cat("second best                        :",
    order(gaps)[2L], " (gap",
    format(sort(gaps)[2L], digits = 3), ")\n")

hr("T3. does a `state` column still get read as known states?")
dh2 <- data.frame(ID = factor(d$id), t = d$t, y = d$y, state = d$state)
hid2 <- suppressMessages(MarkovChain$new(data = dh2, n_states = K,
                                         formula = fmat,
                                         initial_state = "estimated"))
obs2 <- Observation$new(data = dh2, n_states = K, dists = list(y = "norm"),
                        par = list(y = list(mean = c(-2, 0, 3),
                                            sd = rep(0.7, K))))
hm2 <- HMM$new(obs = obs2, hid = hid2)
suppressWarnings(hm2$fit(silent = TRUE))
cat("with a `state` column, hmmTMB llk():", format(hm2$llk(), digits = 12),
    "\n")
cat("without it                         :", format(hm_e$llk(), digits = 12),
    "\n")
cat("difference                         :",
    format(hm2$llk() - hm_e$llk(), digits = 6),
    " <- non-zero means the column changed the model, silently\n")
