# Lane arstep: arstep-r2ng.R with the traces of dev/arstep-trace.R.
# Lane arstep: copy of dev/review-ncgp-r2-ng.R on this round's libraries; arg 3 ref 0 skips the reference starts.
# Fix round 1, item 1: the route with no warm start. Default frm() against the best of a grid of
# starts on the same tape, and on rellib-r16. Families and spellings the gate turns away, on curves
# that have a trend and a wiggle, at high information. Args: arm seeds. Seed 56000 + 13 seed + n.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2])); ref <- if (length(args) > 2) as.integer(args[3]) else 1L
rev_libs(arm)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-trace.R")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x))
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson()),
    pois_high = list(d = data.frame(x = x, y = rpois(n, exp(7 + eta))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli()),
    binom1k = list(d = data.frame(x = x, tr = 1000, y = rbinom(n, 1000, plogis(eta - mean(eta)))), fml = y | trials(tr) ~ gp(x), fam = binomial()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log")),
    lognorm = list(d = data.frame(x = x, y = exp(eta + rnorm(n, sd = 0.02))), fml = y ~ gp(x), fam = lognormal()),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n), y = eta + 0.02 * rnorm(n)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()),
    gmv = list(d = data.frame(x = x, y = eta + 0.02 * rnorm(n), y2 = cos(x) + 0.02 * rnorm(n)), fml = bf(y ~ gp(x)) + bf(y2 ~ gp(x)) + set_rescor(FALSE), fam = NULL),
    glog = list(d = data.frame(x = x, y = exp(eta / 3) + 0.02 * rnorm(n)), fml = y ~ gp(x), fam = gaussian(link = "log")),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
}
fitit <- function(m, start = NULL) {
  a <- list(m$fml, m$d)
  if (!is.null(m$fam)) a$family <- m$fam
  if (!is.null(start)) a$start <- start
  do.call(rev_fit, a)
}
cases <- c("pois_low", "pois_high", "bern", "binom1k", "negbin", "gamma", "lognorm", "student", "gsigma", "gmv", "glog", "gar")
tot <- list()
for (seed in seeds) for (cs in cases) for (cn in names(curves)) for (n in c(60L, 120L)) {
  set.seed(56000 + seed * 13 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(cs, x, curves[[cn]], n)
  sd_log <- character(0); t0 <- proc.time()[[3]]; r <- fitit(m); el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat(sprintf("RNG2 %s %s %s n %d seed %d ERROR %s\n", arm, cs, cn, n, seed, substr(conditionMessage(r$fit), 1, 110))); next }
  ll0 <- as.numeric(logLik(r$fit)); th0 <- r$fit$estimates$theta
  bk <- r$fit$frame$re_blocks
  best <- ll0; bth <- th0
  np <- length(th0)
  if (ref) for (rho0 in c(0.05, 0.15, 0.4, 1, 3, 8)) for (sd0 in c(0.3, 2, 10)) {
    th <- th0
    for (b in bk) if (identical(b$covstruct, "gp")) { th[b$theta_idx[1]] <- log(sd0); th[b$theta_idx[-1]] <- log(rho0) }
    r2 <- fitit(m, list(theta = th))
    if (inherits(r2$fit, "error")) next
    v <- as.numeric(logLik(r2$fit))
    if (is.finite(v) && v > best) { best <- v; bth <- r2$fit$estimates$theta }
  }
  gap <- ll0 - best
  k <- paste(cs)
  tot[[k]] <- rbind(tot[[k]], c(gap, r$fit$opt$convergence, length(r$warn)))
  cat(sprintf("RNG2 %s %s %s n %d seed %d default ll %.6f conv %d w %d theta %s | best %.6f theta %s | gap %.2e %.1fs\n", arm, cs, cn, n, seed, ll0, r$fit$opt$convergence, length(r$warn),
              paste(round(th0, 2), collapse = ","), best, paste(round(bth, 2), collapse = ","), gap, el))
  cat(sprintf("RNG2SD %s %s %s n %d seed %d %s\n", arm, cs, cn, n, seed, if (length(sd_log)) paste(sd_log, collapse = "; ") else "none"))
  flush.console()
}
for (k in names(tot)) { a <- tot[[k]]; cat(sprintf("RNG2SUM %s %-10s rows %d short>1e-3 %d silent %d worst %.3g\n", arm, k, nrow(a), sum(a[, 1] < -1e-3), sum(a[, 1] < -1e-3 & a[, 2] == 0 & a[, 3] == 0), min(a[, 1]))) }
cat("DONE\n")
