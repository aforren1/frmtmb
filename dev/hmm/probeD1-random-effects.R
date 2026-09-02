# Probe D: random effects on state-dependent means. Rung 1 + (1 | g),
# Laplace over b with the forward algorithm inside.
#
# Three references, in increasing independence:
#   (a) a hand-rolled MakeADFun(random = "b") over the same forward
#       algorithm  - same Laplace, so it isolates frmtmb's plumbing;
#   (b) adaptive Gauss-Hermite quadrature over the SCALAR per-group b -
#       the group likelihood given b is exact (forward algorithm), so
#       the marginal is a one-dimensional integral computable to
#       machine precision. This measures the Laplace error itself.
#   (c) hmmTMB, if the model can be spelled there.
#
# Run: FRMTMB_LIB=<lib> Rscript dev/hmm/probeD1-random-effects.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
})
source("dev/hmm/hmm-common.R")
source("dev/hmm/hmm-family.R")

set.seed(97)
K <- 2L
N <- 40L
Tg <- 25L
n <- N * Tg

B_true <- matrix(0, K, K);  B_true[1, 2] <- -1.6; B_true[2, 2] <- 1.2
Bx_true <- matrix(0, K, K)
mu_true <- c(0, 3)
sigma_true <- c(0.6, 0.6)
sd_b <- 0.8
delta_true <- c(0.5, 0.5)

b_true <- rnorm(N, 0, sd_b)
dat <- do.call(rbind, lapply(seq_len(N), function(g) {
  s <- integer(Tg)
  s[1] <- sample.int(K, 1, prob = delta_true)
  G <- tpm_tv_num(B_true, Bx_true, 0, K)
  for (t in seq_len(Tg - 1L)) s[t + 1L] <- sample.int(K, 1, prob = G[s[t], ])
  # the random intercept shifts state 1's mean only
  m <- ifelse(s == 1L, mu_true[1] + b_true[g], mu_true[2])
  data.frame(g = g, t = seq_len(Tg), state = s,
             y = rnorm(Tg, m, sigma_true[s]))
}))
dat$gf <- factor(dat$g)
rows_by_g <- split(seq_len(n), dat$g)

cat("== probe D1: HMM + (1 | g) on state 1's mean ==\n")
cat("   N =", N, " T =", Tg, " sd_b true =", sd_b, "\n\n")

## ---- 1. frm() --------------------------------------------------------

fam <- hmm2_family()
form <- bf(y | vint(g, t) ~ 1 + (1 | gf), mu2 ~ 1, sigma1 ~ 1,
           sigma2 ~ 1, tr12 ~ 1, tr22 ~ 1)
t_fit <- system.time(fit <- frm(form + fam, data = dat))[["elapsed"]]
ll_frm <- as.numeric(logLik(fit))
cat(sprintf("1. frm() fit in %.2f s\n", t_fit))
cat("   logLik :", format(ll_frm, digits = 12), "\n")
fx <- fixef(fit)
sdb_frm <- as.data.frame(VarCorr(fit))$sdcor[1]
cat("   mu1", format(fx$mu[[1]], digits = 6),
    " mu2", format(fx$mu2[[1]], digits = 6),
    " sigma", format(exp(c(fx$sigma1[[1]], fx$sigma2[[1]])), digits = 5),
    "\n")
cat("   tr12", format(fx$tr12[[1]], digits = 6),
    " tr22", format(fx$tr22[[1]], digits = 6),
    " sd(b)", format(sdb_frm, digits = 6), "(true", sd_b, ")\n\n")

## ---- 2. hand-rolled MakeADFun(random = "b") --------------------------

nll <- function(p) {
  "c" <- RTMB::ADoverload("c")
  getAll(p)
  sigma <- exp(lsigma)
  lp <- list(dnorm(dat$y, mu[1] + b[dat$g], sigma[1], log = TRUE),
             dnorm(dat$y, mu[2], sigma[2], log = TRUE))
  eta <- list(list(0, tb[1] + 0 * dat$y), list(0, tb[2] + 0 * dat$y))
  lg <- tpm_logs_ad(eta, K)
  d <- softmax0_ad(ldel, K)
  ll <- sum(dnorm(b, 0, exp(lsdb), log = TRUE))
  for (gi in seq_along(rows_by_g)) {
    ll <- ll + fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], log(d), K)
  }
  -ll
}
par0 <- list(mu = c(fx$mu[[1]], fx$mu2[[1]]),
             lsigma = c(fx$sigma1[[1]], fx$sigma2[[1]]),
             tb = c(fx$tr12[[1]], fx$tr22[[1]]),
             ldel = 0, lsdb = log(sdb_frm), b = rep(0, N))
t_hr <- system.time({
  obj <- MakeADFun(nll, par0, random = "b", silent = TRUE)
  op <- nlminb(obj$par, obj$fn, obj$gr,
               control = list(iter.max = 500, eval.max = 800))
})[["elapsed"]]
cat(sprintf("2. hand-rolled MakeADFun(random = 'b') in %.2f s\n", t_hr))
cat("   logLik :", format(-op$objective, digits = 12), "\n")
cat("   |diff| vs frm :", format(abs(ll_frm + op$objective), digits = 3),
    "\n")
cat("   sd(b)  :", format(exp(op$par[["lsdb"]]), digits = 6),
    " (frm ", format(sdb_frm, digits = 6), ")\n")
cat("   mu     :", format(op$par[1:2], digits = 6), "\n\n")

## ---- 3. Gauss-Hermite quadrature: how good is the Laplace? -----------

# marginal per group: integral over b of exp(fwd(b)) * dnorm(b, 0, sd)
ll_quad <- function(pars, nq) {
  q <- gh_nodes(nq)
  mu <- pars$mu; sigma <- pars$sigma; G <- pars$G
  d <- pars$delta; sdb <- pars$sdb
  tot <- 0
  for (gi in seq_along(rows_by_g)) {
    yg <- dat$y[rows_by_g[[gi]]]
    h <- function(bb) {
      lpm <- cbind(stats::dnorm(yg, mu[1] + bb, sigma[1], log = TRUE),
                   stats::dnorm(yg, mu[2], sigma[2], log = TRUE))
      fwd_num(lpm, G, d) + stats::dnorm(bb, 0, sdb, log = TRUE)
    }
    tot <- tot + aghq1(h, q, -8 * sdb, 8 * sdb)
  }
  tot
}

pe <- op$par
pars_hat <- list(mu = pe[1:2], sigma = exp(pe[3:4]),
                 G = tpm_tv_num(rbind(c(0, pe[5]), c(0, pe[6])),
                                matrix(0, K, K), 0, K),
                 delta = softmax0(pe[["ldel"]]), sdb = exp(pe[["lsdb"]]))
cat("3. adaptive Gauss-Hermite marginal at the Laplace optimum\n")
for (nq in c(3L, 5L, 11L, 21L, 41L)) {
  cat(sprintf("   nq = %3d  logLik = %.10f\n", nq, ll_quad(pars_hat, nq)))
}
llq <- ll_quad(pars_hat, 41L)
cat(sprintf("   Laplace  = %.10f\n", ll_frm))
cat(sprintf("   Laplace - quadrature = %+.3e  (%.2e relative)\n\n",
            ll_frm - llq, abs(ll_frm - llq) / abs(llq)))

## ---- 4. refit by quadrature: parameter bias --------------------------

qnll <- function(v) {
  pl <- list(mu = v[1:2], sigma = exp(v[3:4]),
             G = tpm_tv_num(rbind(c(0, v[5]), c(0, v[6])),
                            matrix(0, K, K), 0, K),
             delta = softmax0(v[7]), sdb = exp(v[8]))
  -ll_quad(pl, 15L)
}
v0 <- c(pe[1:4], pe[5:6], pe[["ldel"]], pe[["lsdb"]])
t_q <- system.time(oq <- optim(v0, qnll, method = "BFGS",
                               control = list(reltol = 1e-12,
                                              maxit = 500)))[["elapsed"]]
cat(sprintf("4. adaptive-quadrature ML refit (nq = 15) in %.1f s\n", t_q))
cat("   logLik  :", format(-oq$value, digits = 12), "\n")
cmp <- rbind(laplace = c(mu1 = pe[1], mu2 = pe[2],
                         sigma1 = exp(pe[3]), sigma2 = exp(pe[4]),
                         tr12 = pe[5], tr22 = pe[6],
                         sdb = exp(pe[["lsdb"]])),
             quad = c(oq$par[1:2], exp(oq$par[3:4]), oq$par[5:6],
                      exp(oq$par[8])))
print(round(cmp, 6))
cat("   max |diff| :", format(max(abs(cmp[1, ] - cmp[2, ])), digits = 3),
    "\n")
cat("   max relative:", format(max(abs(cmp[1, ] - cmp[2, ]) /
                                     pmax(abs(cmp[2, ]), 1e-3)),
                               digits = 3), "\n\n")

## ---- 5. hmmTMB cross-check -------------------------------------------

if (requireNamespace("hmmTMB", quietly = TRUE)) {
  ok <- try({
    suppressPackageStartupMessages(library(hmmTMB))
    dh <- data.frame(ID = dat$g, y = dat$y, gf = dat$gf)
    hid <- MarkovChain$new(data = dh, n_states = 2)
    dists <- list(y = "norm")
    obs <- Observation$new(
      data = dh, dists = dists, n_states = 2,
      par = list(y = list(mean = c(0, 3), sd = c(0.6, 0.6))),
      formulas = list(y = list(mean = ~ s(gf, bs = "re"), sd = ~ 1))
    )
    hm <- HMM$new(obs = obs, hid = hid)
    hm$fit(silent = TRUE)
    hm
  }, silent = TRUE)
  cat("5. hmmTMB\n")
  if (inherits(ok, "try-error")) {
    cat("   FAILED:", substr(attr(ok, "condition")$message, 1, 200), "\n")
  } else {
    cat("   hmmTMB llk :", format(-ok$out()$objective, digits = 12), "\n")
    cat("   frm    llk :", format(ll_frm, digits = 12), "\n")
    cat("   NOTE: hmmTMB puts the RE on BOTH states' means, so the\n",
        "        likelihoods are not the same model; see the memo.\n")
    print(ok$par())
  }
} else {
  cat("5. hmmTMB not installed - skipped\n")
}
cat("\n")

saveRDS(list(dat = dat, ll_frm = ll_frm, ll_hr = -op$objective,
             ll_quad = llq, par_hr = pe, sd_b = sd_b),
        "dev/hmm/probeD1.rds")
cat("saved dev/hmm/probeD1.rds\n")
