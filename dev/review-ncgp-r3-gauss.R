# Short check, items 2 and 3: curves neither the lane nor I have used, both builds, row by row. Reference: my closed form,
# a grid, polished cells, frm()'s point. Seed 73000 + 29 seed + n. Args: arm seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
nug <- if (arm == "lane") 1e-12 else 1e-6
fns <- list(
  ramp = function(x) pmin(x, 4) + 0.3 * sin(5 * x), trig2 = function(x) sin(x) * cos(4 * x), gate = function(x) (x > 2 & x < 6) * 2 + 0.1 * x,
  sqrtw = function(x) 2 * sqrt(x) + 0.5 * cos(6 * x), exp2 = function(x) exp(0.4 * x) + sin(3 * x), kink = function(x) abs(x - 5) + 0.2 * sin(10 * x),
  beat = function(x) sin(3 * x) + sin(3.7 * x), lin = function(x) -30 + 8 * x + cos(2 * x))
for (seed in seeds) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(73000 + seed * 29 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  if (inherits(r$fit, "error")) { cat(sprintf("RG3 %s %s %g %d %d ERROR %s\n", arm, fn, sg, n, seed, substr(conditionMessage(r$fit), 1, 100))); next }
  g <- rev_gp_read(r$fit)
  m <- cf_max(cf_setup(y, x, nug = nug), extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 24L)
  fc <- any(grepl("did not report", r$warn, fixed = TRUE))
  cat(sprintf("RG3 %s %s %g %d %d gap %.3e conv %d fc %d w %d sigma %.4g/%.4g\n", arm, fn, sg, n, seed, g$ll - m$ll, g$conv, fc, length(r$warn), g$sigma, m$sigma))
  flush.console()
}
cat("DONE\n")
