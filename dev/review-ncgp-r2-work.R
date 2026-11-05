# Fix round 1, item 1: what would close the route with no warm start. A working-response start, from outside the
# fit: link(y) with a small clamp, the closed form on it (mine, cf_prof, a coarse grid and polished cells), and its
# sd and length scale handed to frm() as start = list(theta =). Same rows as review-ncgp-r2-ng.R (seed 56000 + 13 seed + n)
# and the spellings of review-ncgp-r2-edge2.R. Reference: `best` from review-ncgp-r2-ng-lane's 18 starts, and this arm's own.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x))
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson(), wk = function(y) log(y + 0.5)),
    pois_high = list(d = data.frame(x = x, y = rpois(n, exp(7 + eta))), fml = y ~ gp(x), fam = poisson(), wk = function(y) log(y + 0.5)),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli(), wk = function(y) qlogis(0.25 + 0.5 * y)),
    binom1k = list(d = data.frame(x = x, tr = 1000, y = rbinom(n, 1000, plogis(eta - mean(eta)))), fml = y | trials(tr) ~ gp(x), fam = binomial(), wk = function(y) qlogis((y + 0.5) / 1001)),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial(), wk = function(y) log(y + 0.5)),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log"), wk = log),
    lognorm = list(d = data.frame(x = x, y = exp(eta + rnorm(n, sd = 0.02))), fml = y ~ gp(x), fam = lognormal(), wk = log),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student(), wk = identity),
    gsigma = list(d = data.frame(x = x, z = rnorm(n), y = eta + 0.02 * rnorm(n)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian(), wk = identity),
    glog = list(d = data.frame(x = x, y = exp(eta / 3) + 0.02 * rnorm(n)), fml = y ~ gp(x), fam = gaussian(link = "log"), wk = function(y) log(pmax(y, 1e-3))),
    gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian(), wk = identity))
}
wstart <- function(yw, x) {
  s <- cf_setup(yw, x)
  d <- sqrt(max(s$D2))
  g <- expand.grid(rho = d / 2^(0:10), lr = log(10^c(1, -1, -3, -6, -9)))
  v <- mapply(function(a, b) cf_prof(s, log(a), b)$ll, g$rho, g$lr)
  rk <- g$rho
  o <- order(-v); o <- o[!duplicated(rk[o])][1:4]
  f <- function(p) { z <- cf_prof(s, p[1], p[2])$ll; if (is.finite(z)) -z else 1e300 }
  best <- NULL
  for (i in o) { q <- nlminb(c(log(g$rho[i]), g$lr[i]), f, control = list(rel.tol = 1e-10)); if (is.null(best) || q$objective < best$objective) best <- q }
  cf_prof(s, best$par[1], best$par[2])
}
ref <- readLines(file.path(REV_OUT, "r2-ng-lane.txt"))
tally <- list()
for (cs in c("pois_low", "pois_high", "bern", "binom1k", "negbin", "gamma", "lognorm", "student", "gsigma", "glog", "gar")) for (cn in names(curves)) for (n in c(60L, 120L)) {
  set.seed(56000 + 13 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(cs, x, curves[[cn]], n)
  t0 <- proc.time()[[3]]; b <- wstart(m$wk(m$d$y), x); tw <- proc.time()[[3]] - t0
  a <- list(m$fml, m$d, family = m$fam, start = list(theta = c(log(b$sd), log(b$rho))))
  t0 <- proc.time()[[3]]; r <- do.call(rev_fit, a); el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat("RWK", cs, cn, n, "ERROR", substr(conditionMessage(r$fit), 1, 100), "\n"); next }
  ll <- as.numeric(logLik(r$fit))
  line <- grep(sprintf("^RNG2 lane %s %s n %d seed 1 ", cs, cn, n), ref, value = TRUE)
  bst <- as.numeric(sub(".*\\| best ([-0-9.]+) theta.*", "\\1", line))
  ll0 <- as.numeric(sub(".*default ll ([-0-9.]+) conv.*", "\\1", line))
  gap <- ll - max(bst, ll)
  tally[[cs]] <- rbind(tally[[cs]], c(gap, ll0 - max(bst, ll)))
  cat(sprintf("RWK %s %s n %d ll %.4f shipped %.4f best-of-18 %.4f | gap %.2e, shipped gap %.2e conv %d w %d start %.2fs fit %.2fs\n", cs, cn, n, ll, ll0, bst, gap, ll0 - max(bst, ll), r$fit$opt$convergence, length(r$warn), tw, el))
  flush.console()
}
for (k in names(tally)) cat(sprintf("RWKSUM %-10s rows %d short: working %d (worst %.3g), shipped %d (worst %.3g)\n", k, nrow(tally[[k]]), sum(tally[[k]][, 1] < -1e-3), min(tally[[k]][, 1]), sum(tally[[k]][, 2] < -1e-3), min(tally[[k]][, 2])))
cat("DONE\n")
