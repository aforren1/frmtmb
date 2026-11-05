# F3 on the built rule: per row, nlminb's code, whether the warning fired, and the gap to
# the closed form with more starts than the fit has. The reviewer's rows (seed 61000 + 23
# seed + n, dev/review-ncgp-r2-fc.R) at seed 1 and fresh at seed 2 and 3, and what
# the rule must leave alone: non-gaussian and a caller's theta.
#   Rscript dev/ncgp-r2fc.R <seeds>
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
seeds <- eval(parse(text = args[1]))
rev_libs("lane")
fns <- list(s2 = function(x) sin(2 * x), s4 = function(x) sin(4 * x) + 0.2 * x, p = function(x) x^2 / 10 - x, m = function(x) 2 * sin(3 * x) + x)
tally <- c(rows = 0, conv1 = 0, warned = 0, quiet = 0, quiet_short = 0, loud_right = 0, loud_short = 0, conv0_short = 0)
worst <- 0
for (seed in seeds) for (fn in names(fns)) for (n in c(30L, 70L, 150L)) for (s in c(1e-2, 1e-3, 3e-4, 1e-4)) {
  set.seed(61000 + 23 * seed + n)
  x <- sort(round(runif(n, 0, 8), 5)); y <- fns[[fn]](x) + s * rnorm(n)
  r <- rev_fit(y ~ gp(x), data.frame(x = x, y = y))
  if (inherits(r$fit, "error")) next
  g <- rev_gp_read(r$fit)
  cs <- cf_setup(y, x)
  m <- cf_max(cs, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 24L)
  gap <- g$ll - m$ll
  wn <- any(grepl("did not report convergence", r$warn, fixed = TRUE))
  tally["rows"] <- tally["rows"] + 1
  if (g$conv != 0) {
    tally["conv1"] <- tally["conv1"] + 1
    if (wn) {
      tally["warned"] <- tally["warned"] + 1
      if (gap < -1e-3) tally["loud_short"] <- tally["loud_short"] + 1 else tally["loud_right"] <- tally["loud_right"] + 1
    } else {
      tally["quiet"] <- tally["quiet"] + 1
      if (gap < -1e-3) { tally["quiet_short"] <- tally["quiet_short"] + 1; worst <- min(worst, gap) }
    }
  } else if (gap < -1e-3) tally["conv0_short"] <- tally["conv0_short"] + 1
  cat(sprintf("RFC %s n %d s=%g seed %d conv %d warned %s gap %.2e\n", fn, n, s, seed, g$conv, wn, gap))
  flush.console()
}
print(tally)
cat("worst among quiet", worst, "\n")
# outside a warm-started gaussian exact gp nothing is quieted
set.seed(61777)
x <- sort(round(runif(60, 0, 8), 5))
y <- sin(2 * x) + 3e-4 * rnorm(60)
d <- data.frame(x = x, y = y, k = rpois(60, exp(1 + sin(x))))
f0 <- suppressWarnings(frm(y ~ gp(x), d))
for (a in list(list(y ~ gp(x), d, start = list(theta = f0$estimates$theta, beta = f0$estimates$beta, betad = f0$estimates$betad)),
               list(k ~ gp(x), d, family = poisson()), list(y ~ x, d))) {
  r <- do.call(rev_fit, a)
  cat(sprintf("RFC2 %s marginal stored %s conv %d warned %s\n", deparse(a[[1]]), !is.null(r$fit$obj$env$.frmtmb_gp_marginal),
              r$fit$opt$convergence, any(grepl("did not report convergence", r$warn, fixed = TRUE))))
}
