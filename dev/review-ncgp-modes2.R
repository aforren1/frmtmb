# Item 1: the remedy and the frequency on a wider set. Arms: the lane; the lane
# started from a local maximize of the closed form (cfopt, the cost is reported);
# Non-gaussian: lane default against the best of several start lengthscales.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- commandArgs(TRUE)[1]
rev_libs(arm)
cfopt <- function(y, x, nug) {
  s <- cf_setup(y, x, nug = nug)
  d <- sqrt(max(s$D2))
  g <- expand.grid(rho = d / 4^(0:5), lr = log(10^c(1, -1, -3, -6, -9)))
  v <- mapply(function(a, b) cf_prof(s, log(a), b)$ll, g$rho, g$lr)
  o <- order(-v)[1:3]
  best <- NULL
  f <- function(p) { z <- cf_prof(s, p[1], p[2])$ll; if (is.finite(z)) -z else 1e300 }
  for (i in o) {
    q <- nlminb(c(log(g$rho[i]), g$lr[i]), f, control = list(rel.tol = 1e-12))
    if (is.null(best) || q$objective < best$objective) best <- q
  }
  cf_prof(s, best$par[1], best$par[2])
}
nug <- if (arm == "lane") 1e-12 else 1e-6
fns <- list(mix = function(x) 2 * sin(3 * x) + x, mix2 = function(x) sin(4 * x) + 0.5 * x,
            mix3 = function(x) 3 * cos(2 * x) - 2 * x, mix4 = function(x) sin(5 * x) + x^2 / 10,
            two = function(x) sin(x) + 0.3 * sin(10 * x), sin3 = function(x) sin(3 * x), sin10 = function(x) sin(10 * x))
tot <- 0; short <- 0; silent <- 0
for (seed in 1:4) for (fn in names(fns)) for (sg in c(0.3, 0.05, 1e-2)) for (n in c(40L, 120L)) {
  set.seed(31000 + seed * 7 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  d <- data.frame(y = y, x = x)
  t0 <- Sys.time(); r <- rev_fit(y ~ gp(x), d); t1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time(); b <- cfopt(y, x, nug); tc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  r2 <- rev_fit(y ~ gp(x), d, start = list(theta = c(log(b$sd), log(b$rho)), betad = log(b$sigma), beta = b$beta[1]))
  t2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  g <- rev_gp_read(r$fit); g2 <- rev_gp_read(r2$fit)
  s <- cf_setup(y, x, nug = nug)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd))), c(log(g2$rho), 2 * (log(g2$sigma) - log(g2$sd)))), ngrid = 26L)
  gap <- g$ll - m$ll; gap2 <- g2$ll - m$ll
  tot <- tot + 1; if (gap < -1e-3) short <- short + 1; if (gap < -1e-3 && g$conv == 0 && !length(r$warn)) silent <- silent + 1
  cat(sprintf("RX %s %s s=%g n %d seed %d cf %.4f sigma %.4g | lane gap %.2e conv %d w %d %.2fs | cfopt gap %.2e conv %d w %d cf %.2fs tape %.2fs\n",
              arm, fn, sg, n, seed, m$ll, m$sigma, gap, g$conv, length(r$warn), t1, gap2, g2$conv, length(r2$warn), tc, t2))
  flush.console()
}
cat(sprintf("RXSUM %s total %d short>1e-3 %d of which silent %d\n", arm, tot, short, silent))
cat("== non-gaussian, default against best of starts\n")
for (seed in 1:3) for (fam in c("poisson", "bernoulli", "binomial")) for (n in c(60L, 150L)) {
  set.seed(32000 + seed + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- 1.5 * sin(3 * x) + 0.4 * x - 2
  d <- data.frame(x = x)
  if (fam == "poisson") { d$y <- rpois(n, exp(eta)); fml <- y ~ gp(x); fm <- poisson() }
  if (fam == "bernoulli") { d$y <- rbinom(n, 1, plogis(eta)); fml <- y ~ gp(x); fm <- bernoulli() }
  if (fam == "binomial") { d$y <- rbinom(n, 50, plogis(eta)); d$t <- 50; fml <- y | trials(t) ~ gp(x); fm <- binomial() }
  r <- rev_fit(fml, d, family = fm)
  if (inherits(r$fit, "error")) { cat("RN", fam, "ERROR", conditionMessage(r$fit), "\n"); next }
  g <- rev_gp_read(r$fit)
  bestll <- -Inf; bestrho <- NA
  for (rho0 in c(0.1, 0.3, 1, 3, 10)) for (sd0 in c(0.5, 3)) {
    r2 <- rev_fit(fml, d, family = fm, start = list(theta = log(c(sd0, rho0))))
    if (inherits(r2$fit, "error")) next
    g2 <- rev_gp_read(r2$fit)
    if (g2$ll > bestll) { bestll <- g2$ll; bestrho <- g2$rho; bestsd <- g2$sd }
  }
  cat(sprintf("RN %s %s n %d seed %d default ll %.4f sd %.4g rho %.4g conv %d w %d | best of 10 starts ll %.4f sd %.4g rho %.4g | gap %.3e\n",
              arm, fam, n, seed, g$ll, g$sd, g$rho, g$conv, length(r$warn), bestll, bestsd, bestrho, g$ll - bestll))
  flush.console()
}
cat("DONE\n")
