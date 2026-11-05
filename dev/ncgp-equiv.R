# B2 on the warm-started build: y and x equivariance across designs, trend plus wiggle
# included. For each design and scale: ll corrected for the unit, sigma, sdgp, lscale,
# fitted SEs at new positions, against the same design at scale 1. Gap to the closed
# form (dev/review-ncgp-lib.R). Seeds in the rows.
#   Rscript dev/ncgp-equiv.R <lib or base>
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- commandArgs(TRUE)[1]
rev_libs(arm)
des <- list(
  sin6 = list(seed = 4140, n = 40L, f = function(x) sin(6 * x), to = 6),
  mix = list(seed = 1, n = 80L, f = function(x) 2 * sin(3 * x) + x, to = 10),
  mix2 = list(seed = 2, n = 60L, f = function(x) sin(4 * x) + 0.5 * x, to = 10),
  two = list(seed = 3, n = 80L, f = function(x) sin(x) + 0.3 * sin(10 * x), to = 10),
  sin2 = list(seed = 20261102, n = 60L, f = function(x) sin(2 * x), to = 6)
)
for (dn in names(des)) for (s in c(0.1, 1e-3, 1e-4)) {
  ds <- des[[dn]]
  set.seed(ds$seed)
  x <- round(stats::runif(ds$n, 0, ds$to), 3)
  e <- stats::rnorm(ds$n)
  y <- ds$f(x) + s * e
  run <- function(ys, xs) {
    d <- data.frame(y = y * ys, x = x * xs)
    r <- rev_fit(y ~ gp(x), d)
    f <- r$fit
    g <- rev_gp_read(f)
    nd <- data.frame(x = xs * seq(min(x) + 0.05, max(x) - 0.07, length.out = 15))
    ft <- suppressWarnings(stats::fitted(f, newdata = nd))
    list(ll = g$ll + ds$n * log(ys), sg = g$sigma / ys, sd = g$sd / ys, rho = g$rho / xs,
         fit = ft[, 1L] / ys, se = ft[, 2L] / ys, conv = g$conv, w = length(r$warn))
  }
  b <- run(1, 1)
  s0 <- cf_setup(y, x)
  cf <- cf_max(s0, extra = list(c(log(b$rho), 2 * (log(b$sg) - log(b$sd)))))
  for (sc in list(c(1, 1), c(1000, 1), c(1e-3, 1), c(1e6, 1), c(1, 1000), c(1, 1e-3))) {
    a <- run(sc[1], sc[2])
    cat(sprintf("EQV %-5s s %.0e y*%-5g x*%-5g conv %d w %d gap %9.2e | ll %9.2e sg %8.1e sd %8.1e rho %8.1e fit %8.1e se %8.1e NaN se %d\n",
                dn, s, sc[1], sc[2], a$conv, a$w, a$ll - cf$ll, a$ll - b$ll, a$sg / b$sg - 1, a$sd / b$sd - 1,
                a$rho / b$rho - 1, max(abs(a$fit - b$fit) / b$se), max(abs(a$se / b$se - 1)), sum(!is.finite(a$se))))
  }
}
