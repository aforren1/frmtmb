# Lane setier, punch round 1: test-perf.R's timed fits, base and lane,
# median of 3 after a warm-up (what that test measures).
#   Rscript dev/setier-perfcheck.R <lib or "base">
args <- commandArgs(TRUE)
.libPaths(c(if (args[1] != "base") args[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
make <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
ft <- function(d) {
  f <- function() suppressWarnings(frm(bf(y ~ x) + poisson(), data = d))
  f()
  stats::median(replicate(3, system.time(f())[["elapsed"]]))
}
cat(args[1], "small", ft(make(1000L, 71)), "large", ft(make(100000L, 72)),
    "\n")
