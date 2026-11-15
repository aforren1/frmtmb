# Lane nlstop2, item 5065 at a larger size: y ~ x + (1 | g), 100 groups of
# `per` rows, t3 noise, seeds 1 to 5. The Laplace fit (quadrature =
# FALSE) against the best of 8 jittered Laplace starts, and the default
# (quadrature on the lane) against 2 jittered starts of its own route, with
# each fit's wall time.
# Usage: NLSTOP2_ARM=after Rscript dev/nlstop2-huber-big.R [per] [seeds]
source("dev/nlstop2-common.R")
a <- commandArgs(TRUE)
per <- if (length(a)) as.numeric(a[1]) else 20
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:5
jit <- function(r, n, quad) {
  best <- llv(r)
  set.seed(n * 31 + 7)
  for (i in seq_len(n)) {
    st <- lapply(r$estimates[c("beta", "betad", "theta")], function(v) {
      unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v)))
    })
    best <- max(best, llv(tt(frm(y ~ x + (1 | g), family = huber(), data = d,
                                 start = st, se = FALSE,
                                 quadrature = quad))), na.rm = TRUE)
  }
  best
}
cat("ROWHEAD route seed ll best gap secs nwarn wclass\n")
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:100, each = per)), x = rnorm(100 * per))
  d$y <- 1 + 0.5 * d$x + 0.5 * rnorm(100)[d$g] + rt(100 * per, 3)
  for (quad in c(FALSE, TRUE)) {
    t0 <- proc.time()[["elapsed"]]
    r <- tt(frm(y ~ x + (1 | g), family = huber(), data = d,
                quadrature = quad))
    t1 <- proc.time()[["elapsed"]] - t0
    b <- jit(r, if (quad) 2 else 8, quad)
    cat(sprintf("ROW %s %d %.6f %.6f %.3g %.2f %d %s\n",
                if (quad) "quadrature" else "laplace", seed, llv(r), b,
                b - llv(r), t1, length(attr(r, "w")), wsum(r)))
    flush.console()
  }
}
cat("DONE\n")
