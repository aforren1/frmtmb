# Lane `latent`, item 2.3, arm B: pin hmmTMB's initial distribution to
# uniform and FIXED, which is what frmtmb's `init = "uniform"` is, then
# ask what the two likelihoods still differ by.
#
#   Rscript dev/latent-2p3-hmmtmb-probe4.R

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
dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"

h0 <- MarkovChain$new(data = dh, n_states = K, formula = fmat,
                      initial_state = "estimated")
ld <- h0$delta0(log = TRUE, as_matrix = FALSE)
cat("log delta0 vector, length", length(ld), ":\n")
print(utils::head(ld, 8))
cat("names:", paste(utils::head(names(ld), 8), collapse = " | "), "\n")

fp <- stats::setNames(rep(NA_real_, length(ld)), rownames(ld))
hid <- MarkovChain$new(data = dh, n_states = K, formula = fmat,
                       initial_state = "estimated",
                       fixpar = list(delta0 = fp))
cat("\nfixpar accepted:\n"); print(utils::head(hid$fixpar()$delta0, 6))

## --- both arms, uniform and fixed initial distribution ---------------
form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fam <- hmm(K = K, gaussian(), time = t, group = id, init = "uniform")
fit <- suppressWarnings(frm(form, family = fam, data = d))
ll_frm <- as.numeric(logLik(fit))
b <- unlist(fixef(fit))
sd_re <- sqrt(VarCorr(fit)[[1L]][1L, 1L])
nre <- nlevels(d$id)

obs <- Observation$new(data = dh, n_states = K, dists = list(y = "norm"),
                       par = list(y = list(mean = c(-2, 0, 3),
                                           sd = rep(0.7, K))))
hm <- HMM$new(obs = obs, hid = hid)
suppressWarnings(hm$fit(silent = TRUE))
ll_h <- hm$llk()
lam <- as.numeric(hm$lambda()$hid)
sd_h <- as.numeric(hid$sd_re()[1L, 1L])

cat("\nfrm  logLik:", format(ll_frm, digits = 12), "\n")
cat("hmm  llk()  :", format(ll_h, digits = 12), "\n")
cat("difference  :", format(ll_frm - ll_h, digits = 8), "\n")
cat("hmmTMB lambda:", format(lam, digits = 8),
    "  1/sqrt(lambda):", format(1 / sqrt(lam), digits = 8),
    "  sd_re:", format(sd_h, digits = 8), "\n")

# The candidate reconstruction: hmmTMB writes the random-effect term as
# a PENALTY, -lambda b'b / 2, with no normalizing constant, while frmtmb
# writes the same term as a normal density. If that is the whole
# difference the two marginals differ by exactly the normalizer.
norm_const <- -(nre / 2) * log(2 * pi * sd_h^2)
cat("\n-(k/2) log(2 pi sd^2) with k =", nre, ":",
    format(norm_const, digits = 8), "\n")
cat("residual after that reconstruction :",
    format(ll_frm - ll_h - norm_const, digits = 8), "\n")

cat("\nestimates\n")
hp <- hm$par()
cf <- hid$coeff_fe()
cat("  means  frm:", format(unname(b[paste0("mu", 1:3, ".(Intercept)")]),
                            digits = 8),
    "\n         tmb:", format(as.numeric(hp$obspar["y.mean", , 1]),
                              digits = 8), "\n")
cat("  sds    frm:", format(exp(unname(
      b[paste0("sigma", 1:3, ".(Intercept)")])), digits = 8),
    "\n         tmb:", format(as.numeric(hp$obspar["y.sd", , 1]),
                              digits = 8), "\n")
cat("  tr12   frm:", format(unname(b["tr12.(Intercept)"]), digits = 8),
    "  tmb:", format(cf["S1>S2.(Intercept)", 1L], digits = 8), "\n")
cat("  tr13   frm:", format(unname(b["tr13.(Intercept)"]), digits = 8),
    "  tmb:", format(cf["S1>S3.(Intercept)", 1L], digits = 8), "\n")
cat("  sd(re) frm:", format(sd_re, digits = 8),
    "  tmb:", format(sd_h, digits = 8), "\n")
