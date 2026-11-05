# F2's claim: the warm start on curves and seeds in neither the review's scripts nor
# test-gp-start.R. Reference: the closed form of dev/review-ncgp-lib.R, a 24-grid plus
# polished cells plus frm()'s own point. Seeds in the row, 71000 + 17 seed + n.
#   Rscript dev/ncgp-r2fresh.R <lane|base> <seeds>
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
nug <- if (arm == "lane") 1e-12 else 1e-6
fns <- list(
  bump3 = function(x) exp(-(x - 2)^2 / 0.5) - exp(-(x - 5)^2 / 0.3) + exp(-(x - 8)^2 / 0.2),
  chirp = function(x) sin(0.4 * x^2),
  tri = function(x) abs((x %% 3) - 1.5) * 2 + 0.3 * x,
  step2 = function(x) (x > 3) - 2 * (x > 7) + 0.2 * sin(6 * x),
  am = function(x) (1 + 0.5 * cos(x)) * sin(5 * x),
  logs = function(x) log1p(x) * 2 + 0.2 * sin(8 * x),
  quad = function(x) (x - 5)^2 / 5 + cos(3 * x),
  off = function(x) 200 + 20 * cos(1.3 * x) + 0.5 * x)
tot <- 0; sh <- 0; sil <- 0; worst <- 0; wall <- 0
for (seed in seeds) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(71000 + seed * 17 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  t0 <- proc.time()[[3]]
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  el <- proc.time()[[3]] - t0
  wall <- wall + el
  if (inherits(r$fit, "error")) { cat("RF3", arm, fn, sg, n, seed, "ERROR", substr(conditionMessage(r$fit), 1, 100), "\n"); next }
  g <- rev_gp_read(r$fit)
  s <- cf_setup(y, x, nug = nug)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 24L)
  gap <- g$ll - m$ll
  tot <- tot + 1
  if (gap < -1e-3) { sh <- sh + 1; if (g$conv == 0 && !length(r$warn)) sil <- sil + 1 }
  worst <- min(worst, gap)
  cat(sprintf("RF3 %s %s s=%g n %d seed %d gap %.2e conv %d w %d sigma %.4g/%.4g %.2fs\n", arm, fn, sg, n, seed, gap, g$conv, length(r$warn), g$sigma, m$sigma, el))
  flush.console()
}
cat(sprintf("RF3SUM %s rows %d short>1e-3 %d silent %d worst %.3g wall %.0fs\n", arm, tot, sh, sil, worst, wall))
