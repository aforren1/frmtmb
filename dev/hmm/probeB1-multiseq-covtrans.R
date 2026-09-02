# Probe B1: N sequences, covariate-dependent transition probabilities.
#
# Question: does the recursion still tape when the transition matrix
# varies by row, and does the resulting ML fit match depmixS4 (which
# supports transition covariates) and a hand-rolled numeric forward?
#
# Run: Rscript dev/hmm/probeB1-multiseq-covtrans.R

suppressPackageStartupMessages(library(RTMB))
source("dev/hmm/hmm-common.R")

set.seed(4321)

K <- 2L
N <- 30L
Tg <- 20L
n <- N * Tg

# transition: row i is a multinomial logit with state 1 as reference.
# Only cell (i, 2) is free for K = 2.
B_true <- matrix(0, K, K)
Bx_true <- matrix(0, K, K)
B_true[1, 2] <- -2.0;  Bx_true[1, 2] <- 1.2   # 1 -> 2
B_true[2, 2] <- 1.5;   Bx_true[2, 2] <- -0.8  # 2 -> 2
mu_true <- c(0, 3)
sigma_true <- c(0.7, 0.7)
delta_true <- c(0.6, 0.4)

dat <- do.call(rbind, lapply(seq_len(N), function(g) {
  x <- rnorm(Tg)
  s <- integer(Tg)
  s[1] <- sample.int(K, 1, prob = delta_true)
  for (t in seq_len(Tg - 1L)) {
    G <- tpm_tv_num(B_true, Bx_true, x[t], K)
    s[t + 1L] <- sample.int(K, 1, prob = G[s[t], ])
  }
  data.frame(g = g, t = seq_len(Tg), x = x, state = s,
             y = rnorm(Tg, mu_true[s], sigma_true[s]))
}))
rows_by_g <- split(seq_len(n), dat$g)

cat("== probe B1: N =", N, "sequences of T =", Tg,
    ", covariate transitions ==\n")
cat("   state occupancy:", table(dat$state), "\n\n")

## ---- taped objective -------------------------------------------------

par0 <- list(mu = c(-0.5, 2.5), lsigma = log(c(1, 1)),
             tb = c(-1, 0.5), tbx = c(0, 0), ldel = 0)

nll <- function(p) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  getAll(p)
  sigma <- exp(lsigma)
  lp <- lapply(seq_len(K), function(k) dnorm(dat$y, mu[k], sigma[k],
                                             log = TRUE))
  # eta[[i]][[j]]: linear predictor for cell (i, j); column 1 is the
  # reference and is never read
  eta <- lapply(seq_len(K), function(i) {
    row <- vector("list", K)
    for (j in seq_len(K)) {
      row[[j]] <- if (j == 1L) 0 else tb[i] + tbx[i] * dat$x
    }
    row
  })
  lg <- tpm_logs_ad(eta, K)
  d <- softmax0_ad(ldel, K)
  ll <- 0
  for (gi in seq_along(rows_by_g)) {
    ll <- ll + fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], log(d), K)
  }
  -ll
}

t_build <- system.time(obj <- MakeADFun(nll, par0, silent = TRUE))[["elapsed"]]
t_opt <- system.time(op <- nlminb(obj$par, obj$fn, obj$gr,
                                  control = list(iter.max = 500,
                                                 eval.max = 800)))[["elapsed"]]
cat("1. taped fit\n")
cat(sprintf("   tape %.2f s  opt %.2f s  conv %d  iters %d\n",
            t_build, t_opt, op$convergence, op$iterations))
cat("   logLik :", format(-op$objective, digits = 12), "\n")
pe <- op$par
cat("   mu     :", format(pe[1:2], digits = 6),
    "(true", mu_true, ")\n")
cat("   sigma  :", format(exp(pe[3:4]), digits = 6),
    "(true", sigma_true, ")\n")
cat("   tb     :", format(pe[5:6], digits = 6),
    "(true", B_true[1, 2], B_true[2, 2], ")\n")
cat("   tbx    :", format(pe[7:8], digits = 6),
    "(true", Bx_true[1, 2], Bx_true[2, 2], ")\n")
cat("   delta  :", format(softmax0(pe[9]), digits = 5),
    "(true", delta_true, ")\n\n")

## ---- numeric reference at the taped optimum --------------------------

B_e <- matrix(0, K, K);  B_e[1, 2] <- pe[5];  B_e[2, 2] <- pe[6]
Bx_e <- matrix(0, K, K); Bx_e[1, 2] <- pe[7]; Bx_e[2, 2] <- pe[8]
lpm <- lpmat_gauss(dat$y, pe[1:2], exp(pe[3:4]))
Gof <- function(r) tpm_tv_num(B_e, Bx_e, dat$x[r], K)
ll_ref <- sum(vapply(rows_by_g, function(rw)
  fwd_num_tv(lpm, Gof, rw, softmax0(pe[9])), numeric(1)))
cat("2. numeric forward at the taped optimum\n")
cat("   taped   :", format(-op$objective, digits = 12), "\n")
cat("   numeric :", format(ll_ref, digits = 12), "\n")
cat("   |diff|  :", format(abs(-op$objective - ll_ref), digits = 3), "\n\n")

## ---- gradient check --------------------------------------------------

g <- as.vector(obj$gr(pe))
h <- 1e-5
fd <- vapply(seq_along(pe), function(i) {
  pp <- pe; pp[i] <- pp[i] + h
  pm <- pe; pm[i] <- pm[i] - h
  (obj$fn(pp) - obj$fn(pm)) / (2 * h)
}, numeric(1))
cat("3. max |AD - FD| gradient:", format(max(abs(g - fd)), digits = 3),
    "\n\n")

## ---- depmixS4 --------------------------------------------------------

if (requireNamespace("depmixS4", quietly = TRUE)) {
  suppressPackageStartupMessages(library(depmixS4))
  set.seed(11)
  best <- NULL
  for (rep in 1:6) {
    m <- depmix(y ~ 1, nstates = K, data = dat, family = gaussian(),
                transition = ~ x, ntimes = rep(Tg, N))
    fm <- try(suppressMessages(fit(m, verbose = FALSE,
                                   emcontrol = em.control(
                                     random.start = TRUE, maxit = 2000))),
              silent = TRUE)
    if (!inherits(fm, "try-error") &&
        (is.null(best) || logLik(fm) > logLik(best))) best <- fm
  }
  cat("4. depmixS4 cross-check (transition = ~ x, ntimes = rep(20, 30))\n")
  cat("   depmixS4 logLik :", format(as.numeric(logLik(best)),
                                     digits = 12), "\n")
  cat("   RTMB     logLik :", format(-op$objective, digits = 12), "\n")
  cat("   |diff|          :",
      format(abs(as.numeric(logLik(best)) + op$objective), digits = 3),
      "\n")
  cat("   depmixS4 pars   :\n")
  print(round(getpars(best), 5))
  cat("\n   RTMB in depmixS4 order:\n")
  cat("     init      :", format(softmax0(pe[9]), digits = 5), "\n")
  cat("     tr row 1  : 0", format(pe[5], digits = 5),
      "| 0", format(pe[7], digits = 5), "\n")
  cat("     tr row 2  : 0", format(pe[6], digits = 5),
      "| 0", format(pe[8], digits = 5), "\n")
  cat("     resp      :", format(pe[1], digits = 5),
      format(exp(pe[3]), digits = 5), "|",
      format(pe[2], digits = 5), format(exp(pe[4]), digits = 5), "\n\n")
} else {
  cat("4. depmixS4 not installed - skipped\n\n")
}

saveRDS(list(dat = dat, par = pe, ll = -op$objective, K = K,
             B_true = B_true, Bx_true = Bx_true, mu_true = mu_true,
             sigma_true = sigma_true, delta_true = delta_true),
        "dev/hmm/probeB1.rds")
cat("saved dev/hmm/probeB1.rds\n")
