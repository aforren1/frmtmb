# N5: the cost of the check as a share of a fit, at the cap, and what the cap and
# the tolerance of integrate() do to the total against every row at a tight
# tolerance. Seed 11502. Elapsed time on a loaded machine, so the share is the
# number to read, from three runs of each.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
pr <- function(...) cat(sprintf(...), "\n", sep = "")
ns <- asNamespace("frmtmb")
check <- ns$mi_laplace_note
set.seed(11502)
for (cfg in list(c(200, 60), c(200, 120), c(1000, 300), c(5000, 1500))) {
  n <- cfg[1]; m <- cfg[2]
  x <- rnorm(n); y <- 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2
  d <- data.frame(x = x, y = y); d$z <- 0.5 + 0.3 * y + rnorm(n, 0, 1)
  d$y[sample(n, m)] <- NA
  tall <- numeric(3); tchk <- numeric(3)
  for (k in 1:3) {
    tall[k] <- system.time(a <- suppressWarnings(frm(
      bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) + gaussian(), data = d)))[["elapsed"]]
    tchk[k] <- system.time(a2 <- suppressWarnings(check(a)))[["elapsed"]]
  }
  # sweeps of the objective per value
  env <- a$obj$env; count <- 0L; f0 <- env$f
  env$f <- function(...) { count <<- count + 1L; f0(...) }
  e <- ns$mi_laplace_row_errors(a, a$mi_laplace$rows)
  sweeps <- count / length(a$mi_laplace$rows)
  env$f <- f0
  pr("COST n=%d read=%d: whole fit %.2fs, the check %.2fs (medians of 3), share of the whole %.1f%%; %d rows, %.0f sweeps per row; total %.2f",
     n, m, median(tall), median(tchk), 100 * median(tchk) / median(tall),
     length(a$mi_laplace$rows), sweeps, a$mi_laplace$total)
  if (n <= 1000) {
    allr <- seq_len(m)
    full <- sum(ns$mi_laplace_row_errors(a, allr), na.rm = TRUE)
    pr("   every row at this tolerance: total %.3f; the sample's scaled total %.3f", full, a$mi_laplace$total)
  }
}
