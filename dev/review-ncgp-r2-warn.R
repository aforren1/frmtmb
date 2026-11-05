# Fix round 1, item 4 and 10: the nugget warning and the false-convergence warning on fresh data.
# Truth noise known; closed form at the same model. The warning is false when sigma is within 0.8 to
# 1.25 of the truth. Fresh seeds and x ranges, none of dev/ncgp-nuggetrate.R's. Args: seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
seeds <- eval(parse(text = args[1]))
rev_libs("lane")
fns <- list(s2 = function(x) sin(2 * x), s4 = function(x) sin(4 * x) + 0.2 * x, bump = function(x) exp(-(x - 4)^2), p = function(x) x^2 / 10 - x)
ratios <- c(1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7)
for (seed in seeds) for (fn in names(fns)) for (n in c(25L, 50L, 110L)) for (s in ratios) {
  set.seed(59000 + 19 * seed + n)
  x <- sort(round(runif(n, 0, 8), 5)); e <- rnorm(n); y <- fns[[fn]](x) + s * e
  r <- rev_fit(y ~ gp(x), data.frame(x = x, y = y))
  if (inherits(r$fit, "error")) next
  g <- rev_gp_read(r$fit)
  nw <- sum(grepl("noise has reached the nugget", r$warn, fixed = TRUE))
  fc <- sum(grepl("did not report", r$warn, fixed = TRUE))
  gw <- sum(grepl("gradient", r$warn, fixed = TRUE))
  within <- sum(grepl("at the maximum and their", r$warn, fixed = TRUE))
  short <- sum(grepl("would still gain|point short", r$warn))
  unus <- sum(grepl("curvature at the stop is unusable", r$warn, fixed = TRUE))
  m <- cf_max(cf_setup(y, x), extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 22L)
  cat(sprintf("RW2 %s n %d seed %d s %.0e sigma/s %.3g sigma/sdgp %.2e sd %.3g gap %.2e conv %d msg %s nw %d fc %d within %d short %d unusable %d grad %d nwarn %d\n", fn, n, seed, s, g$sigma / s, g$sigma / g$sd, g$sd, g$ll - m$ll,
              g$conv, gsub(" ", "_", substr(r$fit$opt$message, 1, 18)), nw, fc, within, short, unus, gw, length(r$warn)))
  flush.console()
}
cat("DONE\n")
