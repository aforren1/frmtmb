# Lane nlstop, item 2569: huber() with a group term (lane optrap's
# construction of dev/optrap-p1-huber19.R: 100 rows, 20 groups, t3
# noise), seeds 1 to 20 at the response scales 1 and 1e3. Each default fit
# against the best of 8 started fits (the default's estimates jittered by
# N(0, 0.3), seed seed * 31 + 7). Short = more than 1e-4 below the best;
# silent = no warning. Also a fit without the group term, whose objective
# has no inner Hessian.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-huber.R
source("dev/nlstop-common.R")
cat("ROWHEAD re scale seed ll conv nwarn wclass best gap\n")
for (re in c(TRUE, FALSE)) for (sc in c(1, 1e3)) for (seed in 1:20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  d$y <- (1 + 0.5 * d$x + 0.5 * b + rt(100, 3)) * sc
  fo <- if (re) bf(y ~ x + (1 | g)) else bf(y ~ x)
  r <- tt(frm(fo, family = huber(), data = d))
  a <- llv(r)
  best <- a
  if (!iserr(r)) {
    set.seed(seed * 31 + 7)
    for (i in 1:8) {
      st <- lapply(r$estimates[c("beta", "betad", "theta")[
        c("beta", "betad", "theta") %in% names(r$estimates)]],
        function(v) unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v))))
      v <- llv(tt(frm(fo, family = huber(), data = d, start = st,
                      se = FALSE)))
      if (is.finite(v) && v > best) best <- v
    }
  }
  cat(sprintf("ROW %s %g %d %.6f %s %d %s %.6f %.3g\n", re, sc, seed, a,
              if (iserr(r)) NA else r$opt$convergence,
              length(attr(r, "w")), wsum(r), best, best - a))
  flush.console()
}
cat("DONE\n")
