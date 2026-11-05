# Fix round 1, item 1: fresh curves and seeds, none in test-gp-start.R. Reference: my closed form
# (dev/review-ncgp-lib.R), a grid plus polished cells plus frm()'s point. Args: arm seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
nug <- if (arm == "lane") 1e-12 else 1e-6
fns <- list(
  saw = function(x) (x %% 2) - 1 + 0.2 * x,
  abs3 = function(x) abs(sin(x)) * 3 - 0.5 * x,
  dbl = function(x) exp(-(x - 3)^2) + 2 * exp(-(x - 7)^2 / 0.1),
  gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x,
  poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x),
  sq = function(x) sign(sin(1.5 * x)) + 0.1 * x,
  tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2,
  big = function(x) 50 + 10 * sin(2.5 * x) + 3 * x)
tot <- 0; sh <- 0; sil <- 0; worst <- 0; cv1 <- 0; nugw <- 0
for (seed in seeds) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(55000 + seed * 11 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  t0 <- proc.time()[[3]]
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat(sprintf("RF2 %s %s s=%g n %d seed %d ERROR %s\n", arm, fn, sg, n, seed, substr(conditionMessage(r$fit), 1, 120))); next }
  g <- rev_gp_read(r$fit)
  s <- cf_setup(y, x, nug = nug)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 24L)
  gap <- g$ll - m$ll
  tot <- tot + 1; if (gap < -1e-3) { sh <- sh + 1; if (g$conv == 0 && !length(r$warn)) sil <- sil + 1 }
  worst <- min(worst, gap); cv1 <- cv1 + (g$conv != 0)
  nw <- sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE)); nugw <- nugw + (nw > 0)
  cat(sprintf("RF2 %s %s s=%g n %d seed %d gap %.2e conv %d w %d nugw %d sigma %.4g/%.4g sd %.4g/%.4g %.2fs\n", arm, fn, sg, n, seed, gap, g$conv, length(r$warn), nw, g$sigma, m$sigma, g$sd, m$sd, el))
  flush.console()
}
cat(sprintf("RF2SUM %s rows %d short>1e-3 %d silent %d worst %.3g conv!=0 %d nugget-warned %d\n", arm, tot, sh, sil, worst, cv1, nugw))
cat("DONE\n")
