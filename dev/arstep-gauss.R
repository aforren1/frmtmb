# Lane arstep: the review's 192 gaussian curves (dev/review-ncgp-r3-gauss.R,
# seeds 21:22, seed 73000 + 29 seed + n) on this round's libraries. Both arms
# are 0.79.x, so both read the relative nugget 1e-12. Reference: the review's
# closed form, a grid, polished cells, frm()'s point. Adds ll, the closed
# form's value and the fit's wall time. Args: arm seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
nug <- 1e-12
fns <- list(
  ramp = function(x) pmin(x, 4) + 0.3 * sin(5 * x), trig2 = function(x) sin(x) * cos(4 * x), gate = function(x) (x > 2 & x < 6) * 2 + 0.1 * x,
  sqrtw = function(x) 2 * sqrt(x) + 0.5 * cos(6 * x), exp2 = function(x) exp(0.4 * x) + sin(3 * x), kink = function(x) abs(x - 5) + 0.2 * sin(10 * x),
  beat = function(x) sin(3 * x) + sin(3.7 * x), lin = function(x) -30 + 8 * x + cos(2 * x))
for (seed in seeds) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(73000 + seed * 29 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  t0 <- proc.time()[[3]]
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  el <- proc.time()[[3]] - t0
  if (inherits(r$fit, "error")) { cat(sprintf("AG %s %s %g %d %d ERROR %s\n", arm, fn, sg, n, seed, substr(conditionMessage(r$fit), 1, 100))); next }
  g <- rev_gp_read(r$fit)
  m <- cf_max(cf_setup(y, x, nug = nug), extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 24L)
  fc <- any(grepl("did not report", r$warn, fixed = TRUE))
  cat(sprintf("AG %s %s %g %d %d gap %.3e conv %d fc %d w %d sigma %.4g/%.4g ll %.6f cf %.6f wall %.2f\n", arm, fn, sg, n, seed, g$ll - m$ll, g$conv, fc, length(r$warn), g$sigma, m$sigma, g$ll, m$ll, el))
  flush.console()
}
cat("DONE\n")
