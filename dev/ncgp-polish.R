# The mutant `nopolish` survives on the final build. Does gp_marginal_polish() still
# earn its cost? Each row fit with it and without, in one process: gap between the two
# log-likelihoods, and against the closed form with more starts than the fit has.
#   Rscript dev/ncgp-polish.R
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
pol <- ns$gp_marginal_polish
put <- function(f) {
  unlockBinding("gp_marginal_polish", ns)
  assign("gp_marginal_polish", f, envir = ns)
  lockBinding("gp_marginal_polish", ns)
}
fns <- list(
  exp8 = function(x) exp(0.8 * x), exp5 = function(x) exp(0.5 * x),
  sq = function(x) x^2, cub = function(x) x^3 / 10, lin = function(x) 5 * x,
  mix = function(x) 2 * sin(3 * x) + x, mix4 = function(x) sin(5 * x) + x^2 / 10,
  bump3 = function(x) exp(-(x - 2)^2 / 0.5) - exp(-(x - 5)^2 / 0.3) + exp(-(x - 8)^2 / 0.2),
  logs = function(x) log1p(x) * 2 + 0.2 * sin(8 * x),
  off = function(x) 200 + 20 * cos(1.3 * x) + 0.5 * x)
tally <- c(rows = 0, polish_gains = 0, polish_loses = 0, on_short = 0, off_short = 0)
worst_off <- 0; worst_on <- 0; t_on <- 0; t_off <- 0
for (seed in 1:2) for (fn in names(fns)) for (sg in c(1, 0.3, 1e-2, 1e-3)) for (n in c(40L, 100L)) {
  set.seed(5700 + seed * 11 + n)
  x <- sort(runif(n, 0, 10))
  y <- fns[[fn]](x) + sg * rnorm(n)
  d <- data.frame(x = x, y = y)
  put(pol)
  t0 <- proc.time()[[3]]
  a <- rev_fit(y ~ gp(x), d)
  t_on <- t_on + proc.time()[[3]] - t0
  put(function(s, out) out)
  t0 <- proc.time()[[3]]
  b <- rev_fit(y ~ gp(x), d)
  t_off <- t_off + proc.time()[[3]] - t0
  put(pol)
  if (inherits(a$fit, "error") || inherits(b$fit, "error")) { cat("RPOL", fn, sg, n, seed, "ERROR\n"); next }
  ga <- rev_gp_read(a$fit); gb <- rev_gp_read(b$fit)
  cs <- cf_setup(y, x)
  m <- cf_max(cs, extra = list(c(log(ga$rho), 2 * (log(ga$sigma) - log(ga$sd))), c(log(gb$rho), 2 * (log(gb$sigma) - log(gb$sd)))), ngrid = 24L)
  on <- ga$ll - m$ll; off <- gb$ll - m$ll
  tally["rows"] <- tally["rows"] + 1
  if (ga$ll - gb$ll > 1e-3) tally["polish_gains"] <- tally["polish_gains"] + 1
  if (gb$ll - ga$ll > 1e-3) tally["polish_loses"] <- tally["polish_loses"] + 1
  if (on < -1e-3) tally["on_short"] <- tally["on_short"] + 1
  if (off < -1e-3) tally["off_short"] <- tally["off_short"] + 1
  worst_on <- min(worst_on, on); worst_off <- min(worst_off, off)
  cat(sprintf("RPOL %s s=%g n %d seed %d on %.2e off %.2e | w %d %d conv %d %d\n", fn, sg, n, seed, on, off, length(a$warn), length(b$warn), ga$conv, gb$conv))
  flush.console()
}
print(tally)
cat(sprintf("RPOLSUM worst on %.3g off %.3g wall on %.0fs off %.0fs\n", worst_on, worst_off, t_on, t_off))
