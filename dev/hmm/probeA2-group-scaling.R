# Probe A2: is the cost driven by the total number of time steps, or by
# the number of sequences? The sizing question for the rung-2 design.
#
# Run: Rscript dev/hmm/probeA2-group-scaling.R

suppressPackageStartupMessages(library(RTMB))
source("dev/hmm/hmm-common.R")

set.seed(31337)
K <- 2L
G_true <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
mu_true <- c(0, 3)
sigma_true <- c(0.6, 0.6)

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

make <- function(N, Tg) {
  dat <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- sim_hmm_seq(Tg, G_true, c(0.5, 0.5), mu_true, sigma_true)
    data.frame(g = g, t = seq_len(Tg), y = s$y)
  }))
  rbg <- split(seq_len(nrow(dat)), dat$g)
  f <- function(p) {
    "c" <- RTMB::ADoverload("c")
    getAll(p)
    sigma <- exp(lsigma)
    lp <- list(dnorm(dat$y, mu[1], sigma[1], log = TRUE),
               dnorm(dat$y, mu[2], sigma[2], log = TRUE))
    eta <- list(list(0, tb[1] + 0 * dat$y), list(0, tb[2] + 0 * dat$y))
    lg <- tpm_logs_ad(eta, K)
    d <- softmax0_ad(ldel, K)
    ll <- 0
    for (k in seq_along(rbg)) {
      ll <- ll + fwd_ad_log_tv(lp, lg, rbg[[k]], log(d), K)
    }
    -ll
  }
  list(f = f, n = nrow(dat))
}

par0 <- list(mu = c(-0.5, 2.5), lsigma = log(c(1, 1)), tb = c(-1.5, 1.5),
             ldel = 0)

cat("== probe A2: cost vs sequence count at fixed total length ==\n\n")
cat(sprintf("   %6s %6s %8s %10s %9s %9s %9s\n", "N", "T", "rows",
            "tape_ms", "fn_ms", "gr_ms", "opt_ms"))
cases <- list(c(1, 5000), c(10, 500), c(50, 100), c(250, 20),
              c(1000, 5))
for (cs in cases) {
  m <- make(cs[1], cs[2])
  tb <- bench(o <- MakeADFun(m$f, par0, silent = TRUE), 3L)
  tf <- bench_at(o$fn, o$par, 20L)
  tg <- bench_at(o$gr, o$par, 20L)
  to <- bench(op <- nlminb(o$par, o$fn, o$gr), 1L)
  cat(sprintf("   %6d %6d %8d %10.1f %9.3f %9.3f %9.1f\n", cs[1], cs[2],
              m$n, tb, tf, tg, to))
}
cat("\n")
cat("   Rows are the unit of cost: the per-sequence loop adds only the\n")
cat("   sequence's own initial and terminal reductions.\n")
