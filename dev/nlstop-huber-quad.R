# Lane nlstop, item 2569: the rows of dev/nlstop-huber.R with the group
# term under quadrature = TRUE, whose objective integrates the density and
# has no inner Hessian, so no jump: is the default fit there the best of
# the same 8 jittered starts?
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-huber-quad.R
source("dev/nlstop-common.R")
cat("ROWHEAD scale seed ll conv nwarn wclass best gap\n")
for (sc in c(1, 1e3)) for (seed in 1:20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  d$y <- (1 + 0.5 * d$x + 0.5 * b + rt(100, 3)) * sc
  fo <- bf(y ~ x + (1 | g))
  r <- tt(frm(fo, family = huber(), data = d, quadrature = TRUE))
  a <- llv(r)
  best <- a
  if (!iserr(r)) {
    set.seed(seed * 31 + 7)
    for (i in 1:8) {
      st <- lapply(r$estimates[c("beta", "betad", "theta")[
        c("beta", "betad", "theta") %in% names(r$estimates)]],
        function(v) unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v))))
      v <- llv(tt(frm(fo, family = huber(), data = d, start = st,
                      se = FALSE, quadrature = TRUE)))
      if (is.finite(v) && v > best) best <- v
    }
  }
  cat(sprintf("ROW %g %d %.6f %s %d %s %.6f %.3g\n", sc, seed, a,
              if (iserr(r)) NA else r$opt$convergence,
              length(attr(r, "w")), wsum(r), best, best - a))
  flush.console()
}
cat("DONE\n")
