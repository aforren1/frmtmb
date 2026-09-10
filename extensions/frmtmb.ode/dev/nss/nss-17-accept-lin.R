# The acceptance criterion, through the closed form.
#
# `frm_lincmt(n_ss = 20)` writes out exactly the twenty cycles
# `frm_ode()` truncates at, and `frm_lincmt()` at its default is the
# limit `frm_ode(ss_extrapolate = TRUE)` now reaches. The two paths
# agree to seven significant digits on this design family
# (dev/lincmt-findings.md), so this is the same contrast at seconds per
# fit instead of tens of minutes, and nss-10-accept.R confirms one pair
# through frm_ode() itself.
#
# Design: two-compartment oral, ke 0.15, k12 0.3, k21 0.02, ka 1.0,
# V 10, ii 24, terminal half-life 107.1 h, 30 subjects x 7 samples,
# random effects on lke and lka, data simulated from the EXACT steady
# state. Seeds 101 to 106.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-17-accept-lin.R
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()

KE <- 0.15; K12 <- 0.3; K21 <- 0.02; KA <- 1.0; V <- 10; II <- 24
NS <- 30L
b <- KE + K12 + K21
LZ <- (b - sqrt(b * b - 4 * KE * K21)) / 2
cat("\nterminal half-life", format(log(2) / LZ), "h, lambda_z * ii",
    format(LZ * II), "\n\n")

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
  list(d = d, ev = ev)
}

fit <- function(z, nss) {
  doses <- z$ev
  main <- if (is.finite(nss))
    conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                   k21 = exp(lk21), ka = exp(lka),
                                   V = exp(lV)),
                      times = time, group = id, ncmt = 2,
                      depot = TRUE, events = doses, n_ss = 20L)
  else
    conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                   k21 = exp(lk21), ka = exp(lka),
                                   V = exp(lV)),
                      times = time, group = id, ncmt = 2,
                      depot = TRUE, events = doses)
  bd <- bf(main, lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
           lk21 ~ 1, lV ~ 1, nl = TRUE)
  frm(bd + gaussian(), data = z$d, se = TRUE,
      start = list(beta = c(log(KE), log(KA), log(K12), log(K21),
                            log(V))))
}

cat(sprintf("%-12s %5s %9s %8s %20s %8s %8s\n", "arm", "seed", "lk21",
            "se", "95 percent interval", "covers", "k21 x"))
cov <- c(limit = 0L, truncated = 0L)
for (seed in 101:106) {
  z <- sim(seed)
  for (nm in c("limit", "truncated")) {
    f <- tryCatch(fit(z, if (nm == "limit") Inf else 20L),
                  error = function(e) e)
    if (inherits(f, "condition")) {
      cat(nm, seed, "ERROR:", conditionMessage(f), "\n")
      next
    }
    ci <- confint(f)
    j <- grep("^lk21_", rownames(ci))[1L]
    est <- ci[j, "est"]
    lo <- ci[j, "lwr"]
    hi <- ci[j, "upr"]
    ok <- isTRUE(log(K21) >= lo && log(K21) <= hi)
    cov[[nm]] <- cov[[nm]] + as.integer(ok)
    cat(sprintf("%-12s %5d %9.4f %8.4f  [%8.4f, %8.4f] %8s %8.2f\n",
                nm, seed, est, (hi - lo) / (2 * 1.96), lo, hi, ok,
                exp(est) / K21))
  }
}
cat("\ntruth log(k21) =", format(log(K21)), "\n")
cat("covered:", paste(names(cov), cov, sep = " ", collapse = ", "),
    "of 6\n")
