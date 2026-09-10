# The plan's acceptance criterion for item 1.0d: a two-compartment oral
# fit at a 107 h terminal half-life recovers k21 inside its own interval
# at the SHIPPED DEFAULT.
#
# The data are simulated from the EXACT steady state, through the closed
# form, so the only thing separating a fit from the truth is the run-in.
# Both arms are frm_ode() at n_ss = 20, the shipped default, so they
# cost the same 21 solves per group per evaluation and differ only in
# whether the geometric tail is summed or dropped.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-10-accept.R
# Seeds: 101, 102, 103. Design: ke 0.15, k12 0.3, k21 0.02, ka 1.0,
# V 10, ii 24, 30 subjects x 7 samples, terminal half-life 107.1 h.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()

KE <- 0.15; K12 <- 0.3; K21 <- 0.02; KA <- 1.0; V <- 10; II <- 24
NS <- 30L
b <- KE + K12 + K21
LZ <- (b - sqrt(b * b - 4 * KE * K21)) / 2
cat("\nterminal half-life", format(log(2) / LZ), "h, lambda_z * ii",
    format(LZ * II), ", r =", format(exp(-LZ * II)), "\n")

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}

sim <- function(seed, sd_obs = 0.15) {
  set.seed(seed)
  tt <- II * c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1)
  d <- data.frame(id = factor(rep(seq_len(NS), each = length(tt))),
                  time = rep(tt, NS))
  i <- as.integer(d$id)
  lke <- log(KE) + rnorm(NS, 0, 0.2)
  lka <- log(KA) + rnorm(NS, 0, 0.3)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = II,
                   addl = 0L, ss = TRUE)
  mu <- numeric(nrow(d))
  for (j in seq_len(NS)) {
    k <- which(i == j)
    mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = K12,
                                     k21 = K21, ka = exp(lka[[j]]),
                                     V = V),
                        times = d$time[k], ncmt = 2, depot = TRUE,
                        events = ev)
  }
  d$conc <- mu + rnorm(nrow(d), 0, sd_obs)
  d
}

# The literal has to be spliced into the formula: a scalar in the
# calling frame is looked up as a data column and fails the length check
# against the data frame.
fit_ode <- function(d, ext) {
  doses <- data.frame(time = 0, state = 1L, value = 100, ii = II,
                      ss = TRUE)
  main <- eval(bquote(
    conc ~ frm_ode(two_oral, init = list(0, 0, 0),
                   times = time, group = id,
                   parms = list(exp(lke), exp(lk12), exp(lk21),
                                exp(lka)),
                   events = doses, output = 2L, n_ss = 20L,
                   ss_tol = 1, ss_extrapolate = .(ext)) / exp(lV)))
  bd <- bf(main,
           lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
           lk21 ~ 1, lV ~ 1, nl = TRUE)
  frm(bd + gaussian(), data = d, se = TRUE,
      start = list(beta = c(log(KE), log(KA), log(K12), log(K21),
                            log(V))))
}

report <- function(tag, f, seed) {
  ci <- confint(f)
  nm <- rownames(ci)
  j <- grep("^lk21_", nm)[1L]
  est <- ci[j, "est"]
  lo <- ci[j, "lwr"]
  hi <- ci[j, "upr"]
  se <- (hi - lo) / (2 * 1.96)
  cat(sprintf("%-14s seed %d  lk21 %8.4f  se %7.4f  95%% [%8.4f,%8.4f]",
              tag, seed, est, se, lo, hi))
  cat(sprintf("  truth %8.4f  covered %s  k21 factor %6.2f\n",
              log(K21), log(K21) >= lo && log(K21) <= hi,
              exp(est) / K21))
  cat(sprintf("               all: %s\n",
              paste(sprintf("%s=%.4f", nm, ci[, "est"]),
                    collapse = " ")))
  cat(sprintf("               logLik %.4f\n", as.numeric(logLik(f))))
  flush(stdout())
  invisible(c(est = est, se = se, cov = log(K21) >= lo &&
                log(K21) <= hi))
}

res <- list()
for (seed in c(101L)) {
  d <- sim(seed)
  for (ext in c(TRUE, FALSE)) {
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(fit_ode(d, ext), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    tag <- if (ext) "extrapolated" else "truncated"
    if (inherits(f, "condition")) {
      cat(tag, "seed", seed, "ERROR:", conditionMessage(f), "\n")
      next
    }
    cat(sprintf("[%.0f s] ", el))
    res[[paste(tag, seed)]] <- report(tag, f, seed)
    saveRDS(res, paste0("C:/Users/adf44/source/r/frmtmb-wt-nss/",
                        "extensions/frmtmb.ode/dev/nss/nss-10-accept.rds"))
  }
}
cat("\ndone\n")
