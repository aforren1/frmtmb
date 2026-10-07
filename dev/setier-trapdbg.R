# Lane setier, punch round 2: the x1e3 trap fits (seeds 21, 23) of
# dev/setier-rev2-trap.R: the objective along the group log sd.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (s in c(21, 23)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * 1e3
  f <- suppressMessages(frm(y ~ x + (1 | g), data = d))
  p <- f$opt$par
  j <- match("theta_1", ns$outer_par_names(f))
  cat("seed", s, "par", signif(p, 5), "gain_up", ns$se_sd_gain_up(f, 1),
      "\n")
  f0 <- f$obj$fn(p)
  for (v in c(p[j] - 2, p[j] + 1, p[j] + 2, log(c(30, 50, 84, 100, 300)))) {
    q <- p
    q[j] <- v
    cat("  theta", signif(v, 4), "dnll", signif(f$obj$fn(q) - f0, 4), "\n")
  }
}
