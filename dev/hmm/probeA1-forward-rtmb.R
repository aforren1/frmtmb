# Probe A1: bare RTMB forward algorithm, 2-state gaussian HMM.
#
# Question: does the forward recursion tape, does it give the right
# likelihood, and which of the two standard formulations (log-space via
# logspace_add, or Zucchini scaling) tapes cleaner and faster?
#
# Run: Rscript dev/hmm/probeA1-forward-rtmb.R

suppressPackageStartupMessages({
  library(RTMB)
})
source("dev/hmm/hmm-common.R")

set.seed(20260902)

K <- 2L
Gamma_true <- matrix(c(0.90, 0.10,
                       0.20, 0.80), 2, 2, byrow = TRUE)
delta_true <- stat_dist(Gamma_true)
mu_true <- c(-1, 2)
sigma_true <- c(0.5, 0.8)

sim <- sim_hmm_seq(200L, Gamma_true, delta_true, mu_true, sigma_true)
y <- sim$y

cat("== probe A1: 2-state gaussian HMM, T =", length(y), "==\n\n")

## ---- 1. numeric reference agrees with itself -------------------------

lpm <- lpmat_gauss(y, mu_true, sigma_true)
ll_scale <- fwd_num(lpm, Gamma_true, delta_true)
ll_log <- fwd_num_log(lpm, Gamma_true, delta_true)
cat("1. numeric forward at the true parameters\n")
cat("   scaled     :", format(ll_scale, digits = 12), "\n")
cat("   log-space  :", format(ll_log, digits = 12), "\n")
cat("   |diff|     :", format(abs(ll_scale - ll_log), digits = 3), "\n\n")

## ---- 2. the two taped versions ---------------------------------------

# parameter vector: mu (K), log sigma (K), transition logits K*(K-1),
# initial-distribution logits (K-1)
par0 <- list(mu = c(-0.5, 1.5), lsigma = log(c(1, 1)),
             lgam = c(-2, -1), ldel = 0)

mk_nll <- function(y, K, kind) {
  Tlen <- length(y)
  function(p) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    getAll(p)
    sigma <- exp(lsigma)
    lp <- lapply(seq_len(K), function(k) dnorm(y, mu[k], sigma[k],
                                               log = TRUE))
    tr <- tpm_rows_ad(lgam, K)
    d <- softmax0_ad(ldel, K)
    if (kind == "log") {
      -fwd_ad_log(lp, tr$lG, log(d), Tlen, K)
    } else {
      -fwd_ad_scale(lp, tr$G, d, Tlen, K)
    }
  }
}

# system.time / proc.time tick at ~15 ms on Windows and one fn() call
# here is microseconds. Run batches of `reps` until at least `min_s`
# has elapsed, then divide; Sys.time() has sub-ms resolution.
bench <- function(expr, reps, min_s = 0.3) {
  expr <- substitute(expr)
  pe <- parent.frame()
  n <- 0L
  t0 <- Sys.time()
  repeat {
    for (i in seq_len(reps)) eval(expr, pe)
    n <- n + reps
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (el >= min_s) break
  }
  1000 * el / n
}

# TMB caches the last evaluation, so obj$fn(obj$par) in a loop times the
# cache, not the tape. Every timing below moves the parameter each call.
bench_at <- function(f, p0, reps, min_s = 0.3) {
  n <- 0L
  t0 <- Sys.time()
  repeat {
    for (i in seq_len(reps)) f(p0 + 1e-6 * ((n + i) %% 11L))
    n <- n + reps
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (el >= min_s) break
  }
  1000 * el / n
}

fits <- list()
for (kind in c("log", "scale")) {
  f <- mk_nll(y, K, kind)
  t_build <- bench(obj <- MakeADFun(f, par0, silent = TRUE), 20L)
  t_fn <- bench_at(obj$fn, obj$par, 200L)
  t_gr <- bench_at(obj$gr, obj$par, 200L)
  t_opt <- bench(op <- nlminb(obj$par, obj$fn, obj$gr), 20L)
  fits[[kind]] <- list(obj = obj, op = op, t_build = t_build,
                       t_fn = t_fn, t_gr = t_gr, t_opt = t_opt)
}

cat("2. taped forward recursions (T = 200), ms per call\n")
for (kind in names(fits)) {
  fi <- fits[[kind]]
  cat(sprintf("   %-6s nll0 %.10f  logLik %.8f  conv %d  iters %d\n",
              kind, fi$obj$fn(unlist(par0)), -fi$op$objective,
              fi$op$convergence, fi$op$iterations))
  cat(sprintf("          tape %.2f ms  fn %.3f ms  gr %.3f ms  opt %.1f ms\n",
              fi$t_build, fi$t_fn, fi$t_gr, fi$t_opt))
}
cat(sprintf("   logLik difference between formulations: %.3e\n\n",
            abs(fits$log$op$objective - fits$scale$op$objective)))

## ---- 3. taped value equals the numeric forward -----------------------

unpack <- function(p, K) {
  mu <- p[1:K]
  sigma <- exp(p[K + 1:K])
  lg <- matrix(p[2 * K + seq_len(K * (K - 1))], K, K - 1, byrow = TRUE)
  G <- tpm_from_logits(lg, K)
  d <- softmax0(p[2 * K + K * (K - 1) + seq_len(K - 1)])
  list(mu = mu, sigma = sigma, G = G, delta = d)
}

pe <- fits$log$op$par
u <- unpack(pe, K)
ll_ref <- fwd_num(lpmat_gauss(y, u$mu, u$sigma), u$G, u$delta)
cat("3. taped optimum vs numeric forward at the SAME parameters\n")
cat("   taped   :", format(-fits$log$op$objective, digits = 12), "\n")
cat("   numeric :", format(ll_ref, digits = 12), "\n")
cat("   |diff|  :", format(abs(-fits$log$op$objective - ll_ref),
                           digits = 3), "\n\n")

## ---- 4. gradient check vs finite differences -------------------------

for (kind in names(fits)) {
  obj <- fits[[kind]]$obj
  p0 <- fits[[kind]]$op$par
  g <- as.vector(obj$gr(p0))
  h <- 1e-5
  fd <- vapply(seq_along(p0), function(i) {
    pp <- p0; pp[i] <- pp[i] + h
    pm <- p0; pm[i] <- pm[i] - h
    (obj$fn(pp) - obj$fn(pm)) / (2 * h)
  }, numeric(1))
  cat(sprintf("4. %-6s max |AD - FD| gradient: %.3e\n", kind,
              max(abs(g - fd))))
}
cat("\n")

## ---- 5. parameter recovery + estimates -------------------------------

cat("5. estimates (log-space formulation)\n")
cat("   mu     :", format(sort(u$mu), digits = 6),
    " (true", format(sort(mu_true)), ")\n")
cat("   sigma  :", format(u$sigma[order(u$mu)], digits = 6),
    " (true", format(sigma_true[order(mu_true)]), ")\n")
cat("   Gamma  :\n")
ord <- order(u$mu)
print(round(u$G[ord, ord], 5))
cat("   true   :\n")
print(round(Gamma_true[order(mu_true), order(mu_true)], 5))
cat("   delta  :", format(u$delta[ord], digits = 5), "\n\n")

## ---- 6. depmixS4 cross-check -----------------------------------------

if (requireNamespace("depmixS4", quietly = TRUE)) {
  suppressPackageStartupMessages(library(depmixS4))
  dd <- data.frame(y = y)
  set.seed(1)
  best <- NULL
  for (rep in 1:5) {
    m <- depmix(y ~ 1, nstates = K, data = dd, family = gaussian())
    fm <- try(suppressMessages(fit(m, verbose = FALSE, emcontrol =
                                     em.control(random.start = TRUE))),
              silent = TRUE)
    if (!inherits(fm, "try-error")) {
      if (is.null(best) || logLik(fm) > logLik(best)) best <- fm
    }
  }
  cat("6. depmixS4 cross-check\n")
  cat("   depmixS4 logLik :", format(as.numeric(logLik(best)),
                                     digits = 12), "\n")
  cat("   RTMB     logLik :", format(-fits$log$op$objective,
                                     digits = 12), "\n")
  cat("   |diff|          :",
      format(abs(as.numeric(logLik(best)) + fits$log$op$objective),
             digits = 3), "\n")
  pp <- getpars(best)
  cat("   depmixS4 pars   :", format(pp, digits = 5), "\n")
  cat("   (order: initial probs, then per-state transition rows,\n",
      "    then per-state (intercept, sd))\n\n")
} else {
  cat("6. depmixS4 not installed - skipped\n\n")
}

## ---- 7. scaling with T -----------------------------------------------

cat("7. tape build / eval cost vs T (log-space formulation)\n")
cat(sprintf("   %7s %10s %10s %10s %10s %8s\n", "T", "tape_ms", "fn_ms",
            "gr_ms", "opt_ms", "nodes"))
for (Tl in c(200L, 1000L, 5000L, 20000L, 50000L)) {
  s2 <- sim_hmm_seq(Tl, Gamma_true, delta_true, mu_true, sigma_true)
  f <- mk_nll(s2$y, K, "log")
  nb <- if (Tl <= 5000L) 5L else 1L
  tb <- bench(o <- MakeADFun(f, par0, silent = TRUE), nb)
  nr <- max(5L, as.integer(2e5 / Tl))
  tf <- bench_at(o$fn, o$par, nr)
  tg <- bench_at(o$gr, o$par, nr)
  to <- bench(op <- nlminb(o$par, o$fn, o$gr), 1L)
  nn <- tryCatch(length(o$env$ADFun$DomainVec()), error = function(e) NA)
  cat(sprintf("   %7d %10.1f %10.3f %10.3f %10.1f %8s\n", Tl, tb,
              tf, tg, to, format(nn)))
}
cat("\n")

## ---- 8. the elementwise-assignment trap ------------------------------

# The gotcha list says a taped loop that assigns single elements is
# >1000x slow. Measure it on the same recursion: build the new alpha
# vector cell by cell instead of folding vectorized logspace_add.
fwd_ad_log_slow <- function(lp, lGrows, ldela, Tlen, K) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  la <- ldela + vapply_ad(lp, 1L, K)
  for (t in seq_len(Tlen - 1L) + 1L) {
    new <- la
    for (j in seq_len(K)) {
      acc <- lGrows[[1]][j] + la[1]
      for (i in seq_len(K - 1L) + 1L) {
        acc <- RTMB::logspace_add(acc, lGrows[[i]][j] + la[i])
      }
      new[j] <- acc + lp[[j]][t]
    }
    la <- new
  }
  acc <- la[1]
  for (i in seq_len(K - 1L) + 1L) acc <- RTMB::logspace_add(acc, la[i])
  acc
}

f_slow <- function(p) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  getAll(p)
  sigma <- exp(lsigma)
  lp <- lapply(seq_len(K), function(k) dnorm(y, mu[k], sigma[k],
                                             log = TRUE))
  tr <- tpm_rows_ad(lgam, K)
  d <- softmax0_ad(ldel, K)
  -fwd_ad_log_slow(lp, tr$lG, log(d), length(y), K)
}
tb_slow <- bench(o_slow <- MakeADFun(f_slow, par0, silent = TRUE), 10L)
tf_slow <- bench_at(o_slow$fn, o_slow$par, 200L)
tg_slow <- bench_at(o_slow$gr, o_slow$par, 200L)
cat("8. elementwise [<-] variant, T = 200\n")
cat(sprintf("   tape %.2f ms (vs %.2f ms vectorized, %.2fx)\n", tb_slow,
            fits$log$t_build, tb_slow / fits$log$t_build))
cat(sprintf("   fn   %.3f ms (vs %.3f ms, %.2fx)\n", tf_slow,
            fits$log$t_fn, tf_slow / fits$log$t_fn))
cat(sprintf("   gr   %.3f ms (vs %.3f ms, %.2fx)\n", tg_slow,
            fits$log$t_gr, tg_slow / fits$log$t_gr))
cat(sprintf("   same value: %.3e\n\n",
            abs(o_slow$fn(unlist(par0)) - fits$log$obj$fn(unlist(par0)))))

saveRDS(list(y = y, par = fits$log$op$par, ll = -fits$log$op$objective,
             Gamma_true = Gamma_true, mu_true = mu_true,
             sigma_true = sigma_true),
        "dev/hmm/probeA1.rds")
cat("saved dev/hmm/probeA1.rds\n")
