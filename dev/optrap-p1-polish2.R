# Lane optrap, punch round 1 (m10): on the degree-6 polynomials, is the
# shortfall the optimizer's? The fit's own objective at lm()'s maximum
# (coefficients and ML sigma) against its stop, and the Newton step's
# predicted and actual change there.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
deg <- 6
for (s in 1:3) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
  f <- suppressWarnings(frm(fo, data = d))
  l <- lm(fo, data = d)
  pl <- c(coef(l), log(sqrt(mean(residuals(l)^2))))
  p <- f$opt$par
  g <- drop(f$obj$gr(p))
  H <- f$obj$he(p)
  st <- -solve(H, g)
  cat(sprintf("seed %d | fn stop %.6f | fn at lm %.6f | -logLik lm %.6f | Newton predicted %.3g actual %.3g | kappa(H) %.2g\n",
              s, f$obj$fn(p), f$obj$fn(pl), -as.numeric(logLik(l)),
              sum(g * st) / 2, f$obj$fn(p + st) - f$obj$fn(p), kappa(H)))
  cat("  stop - lm coefficients / lm SE:",
      signif((p[1:7] - coef(l)) / sqrt(diag(vcov(l))), 2), "\n")
}
