# Short check, item 1: fresh seeds, row by row, both builds. Reference: the best of 15 starts on the same tape, plus the fit's own.
# Curves new in part. Seed 74000 + 31 seed + n. Args: arm seeds cases(comma)
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2])); cases <- strsplit(args[3], ",")[[1]]
rev_libs(arm)
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x), smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson()),
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli()),
    binom5 = list(d = data.frame(x = x, tr = 5, y = rbinom(n, 5, plogis(eta - mean(eta)))), fml = y | trials(tr) ~ gp(x), fam = binomial()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log")),
    lognorm = list(d = data.frame(x = x, y = exp(eta + rnorm(n, sd = 0.02))), fml = y ~ gp(x), fam = lognormal()),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n), y = eta + 0.02 * rnorm(n)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
}
fitit <- function(m, start = NULL) { a <- list(m$fml, m$d, family = m$fam); if (!is.null(start)) a$start <- start; do.call(rev_fit, a) }
for (seed in seeds) for (cs in cases) for (cn in names(curves)) for (n in c(40L, 80L, 160L)) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(cs, x, curves[[cn]], n)
  r <- fitit(m)
  if (inherits(r$fit, "error")) { cat(sprintf("RN3 %s %s %s %d %d ERROR %s\n", arm, cs, cn, n, seed, substr(conditionMessage(r$fit), 1, 100))); next }
  ll0 <- as.numeric(logLik(r$fit)); th0 <- r$fit$estimates$theta; best <- ll0
  for (rho0 in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
    th <- th0
    for (b in r$fit$frame$re_blocks) if (identical(b$covstruct, "gp")) { th[b$theta_idx[1]] <- log(sd0); th[b$theta_idx[-1]] <- log(rho0) }
    r2 <- fitit(m, list(theta = th))
    if (!inherits(r2$fit, "error")) { v <- as.numeric(logLik(r2$fit)); if (is.finite(v) && v > best) best <- v }
  }
  cat(sprintf("RN3 %s %s %s %d %d gap %.3e conv %d w %d ll %.4f best %.4f\n", arm, cs, cn, n, seed, ll0 - best, r$fit$opt$convergence, length(r$warn), ll0, best))
  flush.console()
}
cat("DONE\n")
