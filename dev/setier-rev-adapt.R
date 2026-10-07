# Reviewer of lane setier: the fits behind four test adaptations that now
# allow "Standard errors are not available". Is each lost parameter flat
# at the estimate (the objective moved by +-1 on its internal scale), and
# what did base report for it?
#   Rscript dev/setier-rev-adapt.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.learn)})
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
show <- function(lab, f, probe) {
  nm <- ns$outer_par_names(f)
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf("\n== %s: code %d | lost %s\n", lab, f$opt$convergence,
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none"))
  p <- f$opt$par
  f0 <- f$obj$fn(p)
  for (n in probe) {
    j <- match(n, nm)
    ch <- vapply(c(-1, -0.1, 0.1, 1), function(h) {
      q <- p; q[j] <- q[j] + h; f$obj$fn(q) - f0
    }, 0)
    cat(sprintf("   %-22s est %10.4g SE %10.4g | nll change at -1,-.1,+.1,+1: %s\n",
                n, p[j], se[j], paste(signif(ch, 3), collapse = " ")))
  }
}
set.seed(3)
dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300),
                 g = factor(rep(1:20, 15)))
dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
f <- suppressMessages(suppressWarnings(
  frm(bf(y ~ mo(inc):z + (1 | g)) + gaussian(), data = dm)))
show("brms-likelihood mo(inc):z + (1 | g)", f,
     grep("^zeta|^theta", ns$outer_par_names(f), value = TRUE))

src <- readLines("tests/testthat/helper-reference.R")
eval(parse(text = src[grep("^sim_ar1_data <- function", src):
                        (grep("^sim_ar1_data <- function", src) + 12L)]))
dd <- sim_ar1_data(seed = 53, rho = 0)
f <- suppressMessages(suppressWarnings(
  frm(bf(y ~ 1 + homdiag(tim + 0 | g)) + gaussian(), data = dd)))
show("covstruct homdiag, rho = 0", f,
     grep("^theta|^sigma", ns$outer_par_names(f), value = TRUE))

set.seed(36)
n_g <- 60; n_i <- 8; n <- n_g * n_i
g <- factor(rep(seq_len(n_g), each = n_i))
u <- rnorm(n_g, 0, 0.8); x <- rnorm(n)
lam <- exp(-0.5 + 0.7 * x + u[as.integer(g)])
tt <- stats::rexp(n, lam); ct <- stats::rexp(n, 0.2)
dc <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
f <- suppressMessages(suppressWarnings(
  frm(bf(time | cens(cens) ~ x + (1 | g)), family = cox(), data = dc)))
lost <- names(ns$sdr_of(f)$se_lost)
show("famgaps cox() frailty", f, if (length(lost)) lost else
  ns$outer_par_names(f)[1:2])

d <- frm_task_design("bandit4arm_restless", n_subject = 6L, n_trial = 50L,
                     seed = 56L)
fam <- bandit4arm2_kalman_filter(subject = id, trial = trial)
d$choice <- frm_task_simulate(
  fam, d, pars = list(tau = 0.15, lambda = 0.98, center = 50, mu0 = 50,
                      sigma0 = 12, sigmaD = 2), seed = 56)[[1L]]$choice
f <- suppressMessages(suppressWarnings(frm(
  bf(choice | payoff(pay1, pay2, pay3, pay4) ~ 1, lambda ~ 1, center ~ 1,
     mu0 ~ 1, sigma0 ~ 1, sigmaD ~ 1), family = fam, data = d)))
show("learn Kalman filter", f, ns$outer_par_names(f))
