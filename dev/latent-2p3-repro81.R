# Lane `latent`, item 2.3: reproduce the probe that motivates
# `hmm_starts()`.
#
# dev/hmm-feasibility.md records a cold start converging 8.1
# log-likelihood units below the global optimum with `convergence == 0`,
# `max|grad| == 3.5e-4` and a positive definite Hessian. The
# construction is dev/hmm/probeD4-re-multimodality.R, which predates the
# shipped `hmm()` family and used the rung-1 `custom_family()`
# prototype. This script runs the SAME data and the SAME model through
# the shipped family, so that the figure the help page cites is a figure
# this package can still produce.
#
#   Rscript dev/latent-2p3-repro81.R
#
# Probe D4's data: seed 2026, K = 2, 25 sequences of 30, stationary
# initial distribution, per-state random intercepts.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})

K <- 2L
set.seed(2026)
N <- 25L
Tg <- 30L
G_true <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
mu_true <- c(0, 3)
sigma_true <- c(0.6, 0.6)
sd_b <- c(0.7, 0.5)
stat_dist <- function(G) {
  A <- rbind(t(diag(nrow(G)) - G), 1)
  drop(qr.solve(A, c(rep(0, nrow(G)), 1)))
}
b_true <- cbind(rnorm(N, 0, sd_b[1]), rnorm(N, 0, sd_b[2]))
dat <- do.call(rbind, lapply(seq_len(N), function(g) {
  s <- integer(Tg)
  s[1] <- sample.int(K, 1, prob = stat_dist(G_true))
  for (t in seq_len(Tg - 1L)) {
    s[t + 1L] <- sample.int(K, 1, prob = G_true[s[t], ])
  }
  data.frame(ID = g, t = seq_len(Tg), state = s,
             y = rnorm(Tg, mu_true[s] + b_true[g, s], sigma_true[s]))
}))
dat$gf <- factor(dat$ID)

cat("rows:", nrow(dat), " sequences:", N, " K:", K, "\n")

form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")

t0 <- Sys.time()
fit <- frm(form, family = fam, data = dat)
cat("cold-start fit seconds:",
    as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")
ll_cold <- as.numeric(logLik(fit))
cat("cold start logLik      :", format(ll_cold, digits = 12), "\n")
d <- frmtmb::diagnose(fit, quiet = TRUE)
cat("  convergence", d$convergence, " max|grad|",
    format(d$max_grad, digits = 3), " pdHess", isTRUE(d$pdHess),
    " nbadse", length(d$bad_se), " nflat", length(d$flat), "\n")
cat("  diagnose() message count:", length(d$messages %||% character()),
    "\n")
print(unlist(fixef(fit)))
print(VarCorr(fit))

# probe D4's global optimum, from the hand-rolled MakeADFun arbiter and
# confirmed by hmmTMB to the last digit
ll_global_probe <- -1087.99646521
ll_cold_probe <- -1096.09575602

cat("\nprobe D4 recorded cold start :",
    format(ll_cold_probe, digits = 12), "\n")
cat("probe D4 recorded global     :",
    format(ll_global_probe, digits = 12), "\n")
cat("probe D4 recorded gap        :",
    format(ll_global_probe - ll_cold_probe, digits = 6), "\n")

## Restart at probe D4's global point, expressed in the shipped
## family's own parameter names.
tpl <- par_template(form, family = fam, data = dat)
cat("\npar_template():\n"); print(lapply(tpl, function(v) round(v, 4)))
st <- tpl
st$beta[["mu1_(Intercept)"]] <- -0.185185
st$beta[["mu2_(Intercept)"]] <- 3.121877
st$betad[["sigma1_(Intercept)"]] <- log(0.610455)
st$betad[["sigma2_(Intercept)"]] <- log(0.602054)
st$betad[["tr12_(Intercept)"]] <- log(0.173438 / 0.826562)
st$betad[["tr22_(Intercept)"]] <- log(0.828240 / 0.171760)
st$theta <- log(c(0.645297, 0.427590))
fit2 <- try(frm(form, family = fam, data = dat, start = st), silent = TRUE)
if (inherits(fit2, "try-error")) {
  cat("restart failed:", attr(fit2, "condition")$message, "\n")
} else {
  ll_warm <- as.numeric(logLik(fit2))
  cat("\nrestarted at the probe's global point: logLik",
      format(ll_warm, digits = 12), "\n")
  cat("gap this run reproduces              :",
      format(ll_warm - ll_cold, digits = 6), "\n")
  d2 <- frmtmb::diagnose(fit2, quiet = TRUE)
  cat("  convergence", d2$convergence, " max|grad|",
      format(d2$max_grad, digits = 3), " pdHess", isTRUE(d2$pdHess), "\n")
  print(unlist(fixef(fit2)))
  print(VarCorr(fit2))
}
saveRDS(list(dat = dat, ll_cold = ll_cold), "dev/latent-2p3-repro81.rds")
