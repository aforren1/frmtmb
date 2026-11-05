# N1: the window's Laplace loss per value at the end of the window, on the
# review's shapes (dev/review-refuse-misc.R, seed 10201), a reader of slope 0.2.
# dev/refuse-f1-m6.R integrates one period on the line and so misses the other
# end of the chart; this runs integrate() over the whole chart, as the check does.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(10201)
n <- 200
for (cfg in list(c(3.14, 20), c(3.05, 2), c(0, 2))) {
  ang <- cfg[1] + rnorm(n, 0, 1 / sqrt(cfg[2]))
  y <- ((ang + pi) %% (2 * pi)) - pi
  d <- data.frame(y = y)
  d$z <- 0.2 * y + rnorm(n, 0, 0.5)
  d$y[seq(5, 195, 10)] <- NA
  fit <- suppressWarnings(frm(mvbf(bf(z ~ mi(y)), bf(y | mi() ~ 1, family = von_mises()),
                                   rescor = FALSE), data = d))
  e <- asNamespace("frmtmb")$mi_laplace_row_errors(fit, seq_along(fit$estimates$miss))
  cat(sprintf("mean %.2f kappa %g: %d values, per value mean %.3f range [%.3f, %.3f]; values %.2f %.2f\n",
              cfg[1], cfg[2], length(e), mean(e), min(e), max(e),
              min(frmtmb:::mi_values(fit, "y")[seq(5, 195, 10)]),
              max(frmtmb:::mi_values(fit, "y")[seq(5, 195, 10)])))
}
