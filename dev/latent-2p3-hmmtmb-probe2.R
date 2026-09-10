# Lane `latent`, item 2.3, arm B: WHICH quantity is the identity?
#
# Probe 1 showed the two packages reach the same estimates to about
# 2e-3 while `hm$llk()` and `logLik(fit)` sit 19.7 apart. An identity is
# worthless until the quantity being compared is named, so this script
# asks hmmTMB for every likelihood it can report and finds the one that
# is the same function frmtmb maximizes.
#
#   Rscript dev/latent-2p3-hmmtmb-probe2.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
  library(hmmTMB)
})
source("dev/latent-hmm-sim.R")

s <- hmm_sim(seed = 20260909L, ns = 20L, tl = 100L)
d <- s$d
form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fam <- hmm(K = 3, gaussian(), time = t, group = id, init = "estimated")
fit <- suppressWarnings(frm(form, family = fam, data = d))
ll_frm <- as.numeric(logLik(fit))
cat("frm logLik (Laplace marginal):", format(ll_frm, digits = 12), "\n")
cat("frm df:", attr(logLik(fit), "df"), "\n\n")

dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
K <- 3L
fmat <- matrix("~1", K, K)
diag(fmat) <- "."
fmat[1L, 2L] <- "~s(ID, bs = \"re\")"
hid <- MarkovChain$new(data = dh, n_states = K, formula = fmat,
                       initial_state = "estimated")
b <- unlist(fixef(fit))
obs <- Observation$new(
  data = dh, n_states = K, dists = list(y = "norm"),
  par = list(y = list(mean = c(-2, 0, 3), sd = rep(0.7, K))))
hm <- HMM$new(obs = obs, hid = hid)
suppressWarnings(hm$fit(silent = TRUE))

cat("hm$llk()                      :",
    format(hm$llk(), digits = 12), "\n")
o <- hm$out()
cat("-hm$out()$objective           :",
    format(-o$objective, digits = 12), "  conv", o$convergence, "\n")
to <- hm$tmb_obj()
cat("-tmb_obj$fn(its own par)      :",
    format(-as.numeric(to$fn(to$par)), digits = 12), "\n")
cat("tmb_obj random names          :",
    paste(unique(names(to$env$par[to$env$random])), collapse = ", "),
    "\n")
cat("hm$edf()                      :", format(hm$edf(), digits = 8), "\n")
cat("hm$lambda()                   :\n"); print(hm$lambda())
cat("AIC_marginal                  :", hm$AIC_marginal(), "\n")
cat("AIC_conditional               :", hm$AIC_conditional(), "\n")

cat("\nfrm minus each hmmTMB quantity:\n")
cat("  vs llk()          :", format(ll_frm - hm$llk(), digits = 8), "\n")
cat("  vs -objective     :", format(ll_frm + o$objective, digits = 8),
    "\n")

cat("\nestimates, side by side\n")
hp <- hm$par()
cat("  means  frm:", format(unname(b[c("mu1.(Intercept)",
                                       "mu2.(Intercept)",
                                       "mu3.(Intercept)")]), digits = 8),
    "\n         tmb:", format(as.numeric(hp$obspar["y.mean", , 1]),
                              digits = 8), "\n")
cat("  sds    frm:", format(exp(unname(b[c("sigma1.(Intercept)",
                                           "sigma2.(Intercept)",
                                           "sigma3.(Intercept)")])),
                            digits = 8),
    "\n         tmb:", format(as.numeric(hp$obspar["y.sd", , 1]),
                              digits = 8), "\n")
cf <- hid$coeff_fe()
cat("  tr12   frm:", format(unname(b["tr12.(Intercept)"]), digits = 8),
    "  tmb:", format(cf["S1>S2.(Intercept)", 1L], digits = 8), "\n")
cat("  tr13   frm:", format(unname(b["tr13.(Intercept)"]), digits = 8),
    "  tmb:", format(cf["S1>S3.(Intercept)", 1L], digits = 8), "\n")
cat("  sd(re) frm:", format(sqrt(VarCorr(fit)[[1L]][1L, 1L]), digits = 8),
    "  tmb:", format(as.numeric(hid$sd_re()[1L, 1L]), digits = 8), "\n")
