# Fix round 1, item 10: false convergence on fits at the maximum. For each fit that warns, is the stop at the
# marginal's own maximum, read from gp_marginal_max() on the same frame? That is the proposed quiet rule. Then how
# far short is a fit the rule would silence? Fresh seeds, and a reference with more starts than the fit has.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
seeds <- eval(parse(text = args[1]))
rev_libs("lane")
ns <- asNamespace("frmtmb")
fns <- list(s2 = function(x) sin(2 * x), s4 = function(x) sin(4 * x) + 0.2 * x, p = function(x) x^2 / 10 - x, m = function(x) 2 * sin(3 * x) + x)
tally <- c(rows = 0, fc = 0, within = 0, quiet = 0, quiet_short = 0, loud_right = 0, loud_short = 0)
worst <- 0
for (seed in seeds) for (fn in names(fns)) for (n in c(30L, 70L, 150L)) for (s in c(1e-2, 1e-3, 3e-4, 1e-4)) {
  set.seed(61000 + 23 * seed + n)
  x <- sort(round(runif(n, 0, 8), 5)); y <- fns[[fn]](x) + s * rnorm(n)
  r <- rev_fit(y ~ gp(x), data.frame(x = x, y = y))
  if (inherits(r$fit, "error")) next
  f <- r$fit; g <- rev_gp_read(f)
  fc <- any(grepl("did not report", r$warn, fixed = TRUE))
  within <- any(grepl("at the maximum and their", r$warn, fixed = TRUE))
  mm <- ns$gp_marginal_max(f$frame)
  atm <- !is.null(mm) && abs(f$opt$objective - mm$value) <= 1e-3
  m <- cf_max(cf_setup(y, x), extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 22L)
  gap <- g$ll - m$ll
  quiet <- fc && within && atm
  tally <- tally + c(1, fc, within, quiet, quiet && gap < -1e-3, fc && !quiet && gap >= -1e-3, fc && !quiet && gap < -1e-3)
  if (quiet) worst <- min(worst, gap)
  cat(sprintf("RFC %s n %d seed %d s %.0e conv %d fc %d within %d at-marginal %d stop-minus-marginal %.2e gap %.2e nwarn %d\n", fn, n, seed, s, g$conv, fc, within, atm, -f$opt$objective + (if (is.null(mm)) NA else mm$value), gap, length(r$warn)))
  flush.console()
}
print(tally); cat("worst gap among quieted", worst, "\n")
cat("DONE\n")
