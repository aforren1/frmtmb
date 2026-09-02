# Probe D2: independent third-party cross-check against hmmTMB, on a
# model both packages express exactly: 2-state gaussian HMM, constant
# transitions, STATIONARY initial distribution, with and without
# per-state random intercepts on the state means.
#
# Also the first exercise of the on-tape stationary solve (probe F item).
#
# NOTE: hmmTMB silently reads a data column named `state` as KNOWN
# states. See probeD3. Every hmmTMB call below is given a data frame
# with that column removed.
#
# Run: FRMTMB_LIB=<lib> Rscript dev/hmm/probeD2-hmmtmb.R

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
dat$gf <- factor(dat$ID)
dat_h <- dat[, c("ID", "t", "y", "gf")]   # no `state` column: see probeD3

cat("== probe D2: stationary init + per-state (1 | g), vs hmmTMB ==\n")
cat("   N =", N, " T =", Tg, "\n\n")

## ---- 1. does the on-tape stationary solve work at all? ---------------

fam_s <- hmm2_family_stat()
p_test <- list(mu = rep(0, n), mu2 = rep(3, n), sigma1 = rep(0.6, n),
               sigma2 = rep(0.6, n), tr12 = rep(-1.7, n),
               tr22 = rep(1.4, n))
chk <- tryCatch(
  check_custom_family(fam_s, y = dat$y, dpars = p_test,
                      aterms = list(vint1 = dat$ID, vint2 = dat$t)),
  error = function(e) conditionMessage(e))
cat("1. check_custom_family on the STATIONARY family: ",
    if (isTRUE(chk)) "PASS" else paste("FAIL -", chk), "\n")

G_test <- tpm_tv_num(rbind(c(0, -1.7), c(0, 1.4)), matrix(0, K, K), 0, K)
ll_stat_ref <- sum(vapply(split(seq_len(n), dat$ID), function(r)
  fwd_num_tv(lpmat_gauss(dat$y, c(0, 3), c(0.6, 0.6)),
             function(rr) G_test, r, stat_dist(G_test)), numeric(1)))
ll_stat_fam <- sum(hmm2_lpdf_stat(dat$y, p_test,
                                  list(vint1 = dat$ID, vint2 = dat$t)))
cat("   RTMB::solve stationary delta vs numeric stationary forward:",
    format(abs(ll_stat_fam - ll_stat_ref), digits = 3), "\n\n")

## ---- 2. frm() fits ---------------------------------------------------

form0 <- bf(y | vint(ID, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
            tr12 ~ 1, tr22 ~ 1)
form1 <- bf(y | vint(ID, t) ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf),
            sigma1 ~ 1, sigma2 ~ 1, tr12 ~ 1, tr22 ~ 1)
fit0 <- frm(form0 + fam_s, data = dat)
t_frm <- system.time(fit1 <- frm(form1 + fam_s, data = dat))[["elapsed"]]
ll0 <- as.numeric(logLik(fit0))
ll1 <- as.numeric(logLik(fit1))
fx0 <- fixef(fit0)
fx1 <- fixef(fit1)
vc <- as.data.frame(VarCorr(fit1))
G0 <- tpm_tv_num(rbind(c(0, fx0$tr12[[1]]), c(0, fx0$tr22[[1]])),
                 matrix(0, K, K), 0, K)
G1 <- tpm_tv_num(rbind(c(0, fx1$tr12[[1]]), c(0, fx1$tr22[[1]])),
                 matrix(0, K, K), 0, K)
cat("2. frm()\n")
cat(sprintf("   fixed effects only : logLik %.8f  df %d\n", ll0,
            attr(logLik(fit0), "df")))
cat(sprintf("   + per-state (1|g)  : logLik %.8f  df %d  (%.2f s)\n",
            ll1, attr(logLik(fit1), "df"), t_frm))
cat("   RE model means :", format(c(fx1$mu[[1]], fx1$mu2[[1]]),
                                  digits = 6), " (true", mu_true, ")\n")
cat("   RE model sds   :", format(exp(c(fx1$sigma1[[1]],
                                        fx1$sigma2[[1]])), digits = 5),
    " (true", sigma_true, ")\n")
cat("   sd(b)          :", format(vc$sdcor, digits = 5),
    " (true", sd_b, ")\n")
cat("   tpm            :", format(as.vector(t(G1)), digits = 5),
    " (true", as.vector(t(G_true)), ")\n")
cat("   stationary delta:", format(stat_dist(G1), digits = 5), "\n\n")

## ---- 3. hmmTMB -------------------------------------------------------

if (!requireNamespace("hmmTMB", quietly = TRUE)) {
  cat("3. hmmTMB not installed - skipped\n")
} else {
  suppressPackageStartupMessages(library(hmmTMB))
  run_hmmtmb <- function(re, par0, tpm0) {
    forms <- if (re) {
      list(y = list(mean = ~ s(gf, bs = "re"), sd = ~ 1))
    } else {
      list(y = list(mean = ~ 1, sd = ~ 1))
    }
    hid <- MarkovChain$new(data = dat_h, n_states = 2,
                           initial_state = "stationary", tpm = tpm0)
    obs <- Observation$new(data = dat_h, dists = list(y = "norm"),
                           n_states = 2, par = list(y = par0),
                           formulas = forms)
    hm <- HMM$new(obs = obs, hid = hid)
    suppressWarnings(hm$fit(silent = TRUE))
    hm
  }

  m0 <- c(fx0$mu[[1]], fx0$mu2[[1]])
  s0 <- exp(c(fx0$sigma1[[1]], fx0$sigma2[[1]]))
  hm0 <- run_hmmtmb(FALSE, list(mean = m0, sd = s0), G0)
  ph <- hm0$obs()$par()[, , 1]
  Gh <- hm0$hid()$tpm()[, , 1]
  cat("3. hmmTMB, fixed effects only, both started at frm's optimum\n")
  cat("   frm    logLik :", format(ll0, digits = 12), "\n")
  cat("   hmmTMB llk()  :", format(hm0$llk(), digits = 12), "\n")
  cat("   |diff|        :", format(abs(hm0$llk() - ll0), digits = 3),
      "\n")
  hpar <- c(ph["y.mean", ], ph["y.sd", ])
  cat("   frm par (mu1 mu2 sd1 sd2):", format(c(m0, s0), digits = 8),
      "\n")
  cat("   hmmTMB par               :", format(hpar, digits = 8), "\n")
  cat("   frm tpm   :", format(as.vector(t(G0)), digits = 8), "\n")
  cat("   hmmTMB tpm:", format(as.vector(t(Gh)), digits = 8), "\n")
  cat("   max |par diff|:",
      format(max(abs(c(m0, s0) - hpar),
                 abs(as.vector(t(G0)) - as.vector(t(Gh)))), digits = 3),
      "\n\n")

  hm1 <- run_hmmtmb(TRUE,
                    list(mean = c(fx1$mu[[1]], fx1$mu2[[1]]),
                         sd = exp(c(fx1$sigma1[[1]], fx1$sigma2[[1]]))),
                    G1)
  cat("4. hmmTMB with s(gf, bs = 're') on the mean (RE model)\n")
  cat("   frm    logLik :", format(ll1, digits = 12), "\n")
  cat("   hmmTMB llk()  :", format(hm1$llk(), digits = 12), "\n")
  cat("   |diff|        :", format(abs(hm1$llk() - ll1), digits = 3),
      "\n")
  cat("   frm  sd(b)  :", format(vc$sdcor, digits = 6), "\n")
  sdre <- try(hm1$sd_re(), silent = TRUE)
  if (!inherits(sdre, "try-error")) {
    cat("   hmmTMB sd_re:", format(as.numeric(unlist(sdre$obs)),
                                   digits = 6), "\n")
  }
  ph1 <- hm1$obs()$par()[, , 1]
  cat("   frm means   :", format(c(fx1$mu[[1]], fx1$mu2[[1]]),
                                 digits = 6), "\n")
  cat("   hmmTMB means:", format(ph1["y.mean", ], digits = 6), "\n")
  cat("   frm sds     :", format(exp(c(fx1$sigma1[[1]],
                                       fx1$sigma2[[1]])), digits = 6),
      "\n")
  cat("   hmmTMB sds  :", format(ph1["y.sd", ], digits = 6), "\n")
  cat("   frm tpm     :", format(as.vector(t(G1)), digits = 6), "\n")
  cat("   hmmTMB tpm  :",
      format(as.vector(t(hm1$hid()$tpm()[, , 1])), digits = 6), "\n")
  cat("   hmmTMB -out()$objective:",
      format(-hm1$out()$objective, digits = 12), "\n")
  cat("   hmmTMB AIC_marginal    :",
      format(hm1$AIC_marginal(), digits = 10), "\n")
  cat("   hmmTMB AIC_conditional :",
      format(hm1$AIC_conditional(), digits = 10), "\n")
  cat("   frm AIC                :", format(AIC(fit1), digits = 10),
      "\n")

  ## is frm sitting on a local optimum? restart it at hmmTMB's point
  st <- list(mu = ph1["y.mean", 1], mu2 = ph1["y.mean", 2],
             sigma1 = log(ph1["y.sd", 1]), sigma2 = log(ph1["y.sd", 2]))
  fit1b <- try(frm(form1 + fam_s, data = dat,
                   start = list(beta = c(st$mu, st$mu2, st$sigma1,
                                         st$sigma2, fx1$tr12[[1]],
                                         fx1$tr22[[1]]))), silent = TRUE)
  if (!inherits(fit1b, "try-error")) {
    cat("   frm restarted at hmmTMB's means: logLik",
        format(as.numeric(logLik(fit1b)), digits = 12), "\n")
  } else {
    cat("   frm restart skipped:",
        substr(attr(fit1b, "condition")$message, 1, 120), "\n")
  }
}
cat("\n")
