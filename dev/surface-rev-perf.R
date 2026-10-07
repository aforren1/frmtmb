# Reviewer: test-perf.R's linear-envelope timing, per arm
source("dev/surface-rev-env.R"); arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
make <- function(n, seed) { set.seed(seed); d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(0.3 + 0.4 * d$x)); d }
ft <- function(d) { f <- function() suppressWarnings(frm(bf(y ~ x) + poisson(), data = d))
  f(); median(replicate(3, system.time(f())[["elapsed"]])) }
ts <- ft(make(1000L, 71)); tl <- ft(make(100000L, 72))
cat(sprintf("%s t_small %.3f t_large %.3f ratio %.1f (bound 100)\n", arm, ts, tl, tl / max(ts, 0.01)))
