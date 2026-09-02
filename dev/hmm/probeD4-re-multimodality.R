# Probe D4: frm() and hmmTMB land on different optima for the same
# random-effect HMM (probe D2 section 4). Which likelihood is right, and
# is the gap an implementation difference or plain multimodality?
#
# A hand-rolled MakeADFun(random = "b") over the same model is the
# arbiter: evaluate ITS marginal at both packages' parameter points.
#
# Run: FRMTMB_LIB=<lib> Rscript dev/hmm/probeD4-re-multimodality.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
source("dev/hmm/hmm-common.R")
source("dev/hmm/hmm-family.R")

K <- 2L
set.seed(2026)
N <- 25L
Tg <- 30L
n <- N * Tg
G_true <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
mu_true <- c(0, 3)
sigma_true <- c(0.6, 0.6)
sd_b <- c(0.7, 0.5)
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
rows_by_g <- split(seq_len(n), dat$ID)
gi <- as.integer(factor(dat$ID))

## the same model, hand-rolled: stationary delta, per-state random
## intercepts, constant transitions
nll <- function(p) {
  "c" <- RTMB::ADoverload("c")
  getAll(p)
  sigma <- exp(lsigma)
  sdb <- exp(lsdb)
  lp <- list(dnorm(dat$y, mu[1] + b1[gi], sigma[1], log = TRUE),
             dnorm(dat$y, mu[2] + b2[gi], sigma[2], log = TRUE))
  eta <- list(list(0, tb[1] + 0 * dat$y), list(0, tb[2] + 0 * dat$y))
  lg <- tpm_logs_ad(eta, K)
  g11 <- exp(lg[[1]][[1]][1]); g12 <- exp(lg[[1]][[2]][1])
  g21 <- exp(lg[[2]][[1]][1]); g22 <- exp(lg[[2]][[2]][1])
  G <- RTMB::matrix(c(g11, g21, g12, g22), K, K)
  A <- RTMB::matrix(c(1, 1, 1, 1), K, K) + diag(K) - G
  d <- as.vector(RTMB::solve(t(A), c(1, 1)))
  ll <- sum(dnorm(b1, 0, sdb[1], log = TRUE)) +
    sum(dnorm(b2, 0, sdb[2], log = TRUE))
  for (k in seq_along(rows_by_g)) {
    ll <- ll + fwd_ad_log_tv(lp, lg, rows_by_g[[k]], log(d), K)
  }
  -ll
}

mk <- function(mu, sigma, G, sdb) {
  list(mu = mu, lsigma = log(sigma),
       tb = c(log(G[1, 2] / G[1, 1]), log(G[2, 2] / G[2, 1])),
       lsdb = log(sdb), b1 = rep(0, N), b2 = rep(0, N))
}

# frm() point (probe D2 section 2) and hmmTMB point (section 4)
p_frm <- mk(c(-0.188984, 2.993285), c(0.612995, 0.599946),
            matrix(c(0.811999, 0.188001, 0.173479, 0.826521), 2, 2,
                   byrow = TRUE), c(0.64559, 0.75477))
p_hmm <- mk(c(0.258342, 3.952769), c(0.610454, 0.602054),
            matrix(c(0.826562, 0.173438, 0.171761, 0.828239), 2, 2,
                   byrow = TRUE), c(0.64530, 0.42759))

obj <- MakeADFun(nll, p_frm, random = c("b1", "b2"), silent = TRUE)
pv <- function(p) unlist(p[setdiff(names(p), c("b1", "b2"))])

cat("== probe D4: whose optimum is it? ==\n\n")
cat("hand-rolled Laplace marginal, evaluated at:\n")
cat("   frm's    reported parameters :",
    format(-obj$fn(pv(p_frm)), digits = 12), "\n")
cat("   hmmTMB's reported parameters :",
    format(-obj$fn(pv(p_hmm)), digits = 12), "\n")
cat("   (both transcribed to 6 digits, and hmmTMB's obs()$par() for a\n")
cat("    smooth model reports a fitted value, not the intercept, so\n")
cat("    neither row is the packages' own objective: the optimizations\n")
cat("    below are what matters)\n\n")

for (nm in c("frm", "hmm")) {
  p0 <- if (nm == "frm") p_frm else p_hmm
  op <- nlminb(pv(p0), obj$fn, obj$gr,
               control = list(iter.max = 800, eval.max = 1200))
  cat(sprintf("optimized from the %s start: logLik %.8f  conv %d\n",
              nm, -op$objective, op$convergence))
  cat("   mu   :", format(op$par[1:2], digits = 6),
      "  sigma:", format(exp(op$par[3:4]), digits = 6), "\n")
  cat("   sd(b):", format(exp(op$par[7:8]), digits = 6),
      "  tpm  :", format(as.vector(t(tpm_tv_num(
        rbind(c(0, op$par[5]), c(0, op$par[6])), matrix(0, K, K), 0, K))),
        digits = 6), "\n")
}
cat("\ntrue values: mu", mu_true, " sigma", sigma_true,
    " sd(b)", sd_b, "\n")
cat("true tpm   :", as.vector(t(G_true)), "\n\n")

## ---- is frm()'s objective the same FUNCTION? -------------------------

dat$gf <- factor(dat$ID)
fam_s <- hmm2_family_stat()
form1 <- bf(y | vint(ID, t) ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf),
            sigma1 ~ 1, sigma2 ~ 1, tr12 ~ 1, tr22 ~ 1)
fit1 <- frm(form1 + fam_s, data = dat)

# frm's obj$par is its STARTING vector: beta, betad(5), theta(2), with
# theta the log-SDs. Evaluate both objectives there.
pf <- fit1$obj$par
v <- c(pf[names(pf) == "beta"], pf[names(pf) == "betad"],
       pf[names(pf) == "theta"])
p_at <- list(mu = c(v[1], v[2]), lsigma = c(v[3], v[4]),
             tb = c(v[5], v[6]), lsdb = c(v[7], v[8]),
             b1 = rep(0, N), b2 = rep(0, N))
cat("both objectives at frm's own cold-start vector:\n")
cat("   frm         :", format(-fit1$obj$fn(pf), digits = 14), "\n")
cat("   hand-rolled :", format(-obj$fn(pv(p_at)), digits = 14), "\n")
cat("   |diff|      :",
    format(abs(fit1$obj$fn(pf) - obj$fn(pv(p_at))), digits = 3), "\n\n")

cat("frm from its default cold start : logLik",
    format(as.numeric(logLik(fit1)), digits = 12), "\n")
d <- try(diagnose(fit1), silent = TRUE)
if (!inherits(d, "try-error")) {
  cat("   convergence code", d$convergence, " max|grad|",
      format(d$max_grad, digits = 3), " pdHess", d$pdHess, "\n")
}

## restart frm at the hand-rolled optimum's fixed effects
fit1c <- try(frm(form1 + fam_s, data = dat,
                 start = list(
                   beta = -0.185185,
                   betad = c(3.121877, log(0.610455), log(0.602054),
                             log(0.173438 / 0.826562),
                             log(0.828240 / 0.171760)))),
             silent = TRUE)
if (inherits(fit1c, "try-error")) {
  cat("frm restart error:",
      substr(attr(fit1c, "condition")$message, 1, 200), "\n")
} else {
  cat("frm restarted at the global optimum : logLik",
      format(as.numeric(logLik(fit1c)), digits = 12), "\n")
  cat("   hand-rolled optimum               :  -1087.99646521\n")
  cat("   hmmTMB optimum                    :  -1087.99646521\n")
}
