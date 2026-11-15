.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (cs in list(c(4, 1e-4), c(5, 1e-4), c(1, 1e-5), c(2, 1e-5))) {
  seed <- cs[1]; s <- cs[2]
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * x0
  y <- eta + rnorm(300)
  d <- data.frame(g = g, x = x0 * s, y = y)
  f <- frm(y ~ 1 + (1 + x | g), data = d,
           control = frmtmb_control(autoscale = FALSE))
  d$x <- x0
  fr <- frm(y ~ 1 + (1 + x | g), data = d)
  cat("seed", seed, "s", s, "short", logLik(fr) - logLik(f), "\n")
  print(f$opt$par)
  print(fr$opt$par)
}
