# Lane nlstop2, item 5065: huber() with (1 | g) under the Laplace
# approximation. The rows of dev/nlstop-huber.R (100 rows, 20 groups, t3
# noise, seeds 1 to 20, response scales 1 and 1e3), each default fit
# against the best of 8 jittered starts. A candidate detector: how near
# the kink k a standardized residual sits at the fit, min | |u_i| - k |
# with u_i = (y_i - mu_i) / sigma at the random effects' modes, where the
# Laplace objective jumps.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-huberkink.R [scales]
source("dev/nlstop2-common.R")
a <- commandArgs(TRUE)
scales <- if (length(a)) eval(parse(text = a[1])) else c(1, 1e3)
k <- 1.345
kink <- function(r, d) {
  mu <- as.numeric(fitted(r))
  s <- exp(unname(r$estimates$betad[1]))
  u <- (d$y - mu) / s
  dd <- abs(abs(u) - k)
  c(min = min(dd), n6 = sum(dd < 1e-6), n4 = sum(dd < 1e-4))
}
cat("ROWHEAD scale seed ll conv nwarn wclass best gap t kmin n6 n4\n")
for (sc in scales) for (seed in 1:20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  d$y <- (1 + 0.5 * d$x + 0.5 * b + rt(100, 3)) * sc
  fo <- bf(y ~ x + (1 | g))
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(fo, family = huber(), data = d))
  t1 <- proc.time()[["elapsed"]] - t0
  a <- llv(r)
  best <- a
  km <- c(NA, NA, NA)
  if (!iserr(r)) {
    km <- kink(r, d)
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
  cat(sprintf("ROW %g %d %.6f %s %d %s %.6f %.3g %.3f %.3g %d %d\n", sc, seed,
              a, if (iserr(r)) NA else r$opt$convergence,
              length(attr(r, "w")), wsum(r), best, best - a, t1,
              km[1], km[2], km[3]))
  flush.console()
}
cat("DONE\n")
