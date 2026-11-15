# Short check of lane slopescale's fix round: cost of the default's
# engagement on the dose01 design of slopescale-short-trigger.R (x 0 or
# 1 + 0.01 z), where autoscale = FALSE already reaches the maximum.
# Lane library only; each time is the minimum of 5 blocks grown past 1 s.
# Usage: Rscript slopescale-short-cost.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
mk <- function(seed) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  trt <- rep(c(TRUE, FALSE), length.out = 300)
  z <- rnorm(300)
  x <- ifelse(trt, 1 + 0.01 * z, 0)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + (0.3 + rnorm(20, 0, 0.5)[g]) * x
  data.frame(g = g, x = x, y = eta + rnorm(300))
}
tm <- function(expr) {
  f <- function() eval(expr)
  n <- 1L
  repeat {
    t <- system.time(for (i in seq_len(n)) suppressMessages(f()))[["elapsed"]]
    if (t > 1) break
    n <- n * 2L
  }
  min(replicate(5, system.time(for (i in seq_len(n))
    suppressMessages(f()))[["elapsed"]])) / n
}
for (seed in 1:2) for (fo in c("y ~ x + (1 + x | g)", "y ~ 1 + (1 + x | g)")) {
  d <- mk(seed)
  fo <- as.formula(fo)
  td <- tm(quote(frm(fo, data = d)))
  tf <- tm(quote(frm(fo, data = d,
                     control = frmtmb_control(autoscale = FALSE))))
  cat(sprintf("seed %d %-22s default %.4f s  FALSE %.4f s  ratio %.2f\n",
              seed, deparse(fo), td, tf, td / tf))
}
