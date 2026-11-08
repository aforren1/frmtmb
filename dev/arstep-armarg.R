# Lane arstep: the AR marginal against the fit's objective, and the fits.
# 1. On each row, gp_marginal_nll() at 20 random outer points against
#    fit$obj$fn() at the same points (gaussian: the Laplace step is exact,
#    so the two are one number up to rounding). Seed 9100.
# 2. The default fit against the review's 15 reference starts.
# Rows: the step rows of the seed scan (dev/arstep-seedscan.R), data as
# dev/review-ncgp-r3-ng.R (seed 74000 + 31 seed + n); an ar(p = 2) and a
# grouped ar() row. Args: arm ref(0/1)
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; ref <- if (length(args) > 1) as.integer(args[2]) else 1L
rev_libs(arm)
ns <- asNamespace("frmtmb")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, step = function(x) 2 * (x > 5) - 1)
rows <- list(c("gar", "step", 160, 1), c("gar", "step", 160, 2), c("gar", "step", 80, 22), c("gar", "step", 80, 23),
             c("gar", "step", 80, 4), c("gar", "step", 40, 2), c("gar", "mix", 80, 1),
             c("gar2", "step", 80, 1), c("ggr", "step", 120, 1))
for (rw in rows) {
  cs <- rw[1]; cn <- rw[2]; n <- as.integer(rw[3]); seed <- as.integer(rw[4])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- curves[[cn]](x)
  m <- switch(cs,
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t)),
    gar2 = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = c(0.4, 0.2)), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t, p = 2)),
    ggr = list(d = data.frame(x = x, t = rep(seq_len(n / 3), 3), g = rep(1:3, each = n / 3), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t, gr = g)))
  t0 <- proc.time()[[3]]
  r <- rev_fit(m$fml, m$d, family = gaussian())
  el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat(sprintf("AM %s %s %d %d ERROR %s\n", cs, cn, n, seed, conditionMessage(r$fit))); next }
  fit <- r$fit
  s <- ns$gp_marginal_setup(fit$frame)
  if (!is.null(s)) {
    nll <- ns$gp_marginal_nll(s)
    set.seed(9100)
    obj <- fit$obj
    p0 <- fit$opt$par
    d <- numeric(20)
    for (k in 1:20) {
      p <- p0 + rnorm(length(p0), sd = 0.3)
      pl <- obj$env$parList(p)
      # ML: the tape holds beta in natural coordinates (no whitening)
      stopifnot(is.null(fit$beta_whiten))
      pv <- list(beta = pl$beta, betad = pl$betad, theta = pl$theta, thetaac = pl$thetaac)
      a <- nll(pv)
      b <- obj$fn(p)
      d[k] <- (a - b) / max(1, abs(b))
    }
    cat(sprintf("AM %s %s %d %d marginal-vs-objective max rel diff %.3e over 20 points\n", cs, cn, n, seed, max(abs(d))))
  } else cat(sprintf("AM %s %s %d %d marginal gate refuses\n", cs, cn, n, seed))
  ll0 <- as.numeric(logLik(fit)); best <- ll0
  if (ref) for (rho0 in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
    th <- fit$estimates$theta; th[1] <- log(sd0); th[2] <- log(rho0)
    r2 <- rev_fit(m$fml, m$d, family = gaussian(), start = list(theta = th))
    if (!inherits(r2$fit, "error")) { v <- as.numeric(logLik(r2$fit)); if (is.finite(v) && v > best) best <- v }
  }
  e <- fit$estimates
  cat(sprintf("AM %s %s %d %d %s ll %.6f best %.6f gap %.3e conv %d w %d wall %.2f sdgp %.4g lscale %.4g ar %s\n", cs, cn, n, seed, arm, ll0, best, ll0 - best,
              fit$opt$convergence, length(r$warn), el, exp(e$theta[1]), exp(e$theta[2]), paste(signif(e$thetaac, 4), collapse = ",")))
  if (length(r$warn)) cat("AMW", substr(r$warn, 1, 120), "\n")
  flush.console()
}
cat("DONE\n")
