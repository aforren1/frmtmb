# Stage 4: random effects in state dpars; hmmTMB agreement.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")

## probe D2's data, verbatim
K <- 2L
set.seed(2026)
N <- 25L; Tg <- 30L; n <- N * Tg
G_true <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
mu_true <- c(0, 3); sigma_true <- c(0.6, 0.6); sd_b <- c(0.7, 0.5)
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

fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit0 <- frm(bf(y ~ 1), family = fam, data = dat)
cat("fixed-effects logLik:", sprintf("%.8f", as.numeric(logLik(fit0))),
    "\nprobe D2 reference  : -1216.40337453 (hmmTMB 8.7e-10)\n")
e <- unlist(fixef(fit0))
cat("mu   :", round(c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]]), 8),
    "\nsd   :", round(exp(c(e[["sigma1.(Intercept)"]],
                           e[["sigma2.(Intercept)"]])), 8), "\n")
lg <- c(e[["tr12.(Intercept)"]], e[["tr22.(Intercept)"]])
G <- rbind(c(1, exp(lg[1])) / (1 + exp(lg[1])),
           c(1, exp(lg[2])) / (1 + exp(lg[2])))
cat("tpm  :", round(as.vector(t(G)), 8), "\n")
cat("probe: -0.18119852 3.07922909 | 0.80618908 0.79517108 |",
    "0.83056236 0.16943764 0.15862241 0.84137759\n\n")

## --- hmmTMB, started at frm's optimum --------------------------------
if (requireNamespace("hmmTMB", quietly = TRUE)) {
  dat_h <- dat[, c("ID", "t", "y")]   # NEVER pass `state`: probeD3
  hid <- hmmTMB::MarkovChain$new(data = dat_h, n_states = 2,
                                 initial_state = "stationary")
  obs <- hmmTMB::Observation$new(
    data = dat_h, n_states = 2, dists = list(y = "norm"),
    par = list(y = list(mean = c(e[["mu1.(Intercept)"]],
                                 e[["mu2.(Intercept)"]]),
                        sd = exp(c(e[["sigma1.(Intercept)"]],
                                   e[["sigma2.(Intercept)"]])))))
  hid$update_tpm(G)
  hm <- hmmTMB::HMM$new(obs = obs, hid = hid)
  hm$fit(silent = TRUE)
  cat("hmmTMB logLik:", sprintf("%.8f", hm$llk()),
      " diff:", abs(hm$llk() - as.numeric(logLik(fit0))), "\n")
  hp <- hm$par()
  print(hp$obspar)
  print(hp$tpm)
}

## --- random effects ---------------------------------------------------
t1 <- system.time(
  fit1 <- frm(bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf)),
              family = fam, data = dat))[["elapsed"]]
cat("\nRE logLik:", sprintf("%.8f", as.numeric(logLik(fit1))),
    " df", attr(logLik(fit1), "df"), " time", round(t1, 2), "s\n")
print(VarCorr(fit1))
cat("\nprobe D4 global optimum: -1087.99646521 with sigma ",
    "0.610455/0.602054 and sd(b) 0.645297/0.427590\n")
e1 <- unlist(fixef(fit1))
cat("sigma:", round(exp(c(e1[["sigma1.(Intercept)"]],
                          e1[["sigma2.(Intercept)"]])), 6), "\n")

## --- quadrature refusal ----------------------------------------------
r <- try(frm(bf(y ~ 1 + (1 | gf)), family = fam, data = dat,
             quadrature = TRUE), silent = TRUE)
cat("\nquadrature refusal:", conditionMessage(attr(r, "condition")), "\n")
r <- try(frm(bf(y ~ 1), family = fam, data = dat, REML = TRUE),
         silent = TRUE)
cat("\nREML refusal:", conditionMessage(attr(r, "condition")), "\n")
