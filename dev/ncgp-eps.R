# The step of frm_curve_deriv() on an exact gp() at relative 1e-12: errors against
# the kernel's own derivative, which is exact. Seeds 5 and 17. Prints, per step in
# units of the grid's range, the worst relative error of the derivative's value and
# of its standard error against the largest step's.
#   Rscript dev/ncgp-eps.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
des <- list(
  a = local({set.seed(5); x <- sort(stats::runif(50, 0, 4)); data.frame(x = x, y = sin(1.3 * x) + stats::rnorm(50, 0, 0.15))}),
  b = local({set.seed(17); x <- round(stats::runif(120, 0, 10), 1); data.frame(x = x, y = 2 + sin(x) + stats::rnorm(120, 0, 0.5))}),
  c = local({set.seed(17); x <- stats::runif(200, 0, 10); data.frame(x = x, y = sin(2 * x) + stats::rnorm(200, 0, 0.01))}))
for (dn in names(des)) {
  d <- des[[dn]]
  fg <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
  th <- fg$estimates$theta
  s <- exp(th[1]); l <- exp(th[2])
  pos <- fg$frame$linpreds[["y.mu"]]$gps[[1]]$positions[, 1]
  K <- exp(-outer(pos, pos, "-")^2 / (2 * l^2)) + diag(1e-12, length(pos))
  R <- chol(K)
  z <- fg$estimates$b
  xs <- seq(min(d$x) + 0.1, max(d$x) + 0.5, length.out = 15)
  nd <- data.frame(x = xs)
  r <- outer(xs, pos, "-"); k0 <- exp(-r^2 / (2 * l^2))
  dk <- list(-r / l^2 * k0, (r^2 / l^4 - 1 / l^2) * k0)
  rg <- diff(range(xs))
  for (o in 1:2) {
    ex <- s * drop(t(backsolve(R, t(dk[[o]]), transpose = TRUE)) %*% z)
    ref <- frm_curve_deriv(fg, var = "x", order = o, newdata = nd, simultaneous = FALSE, eps = rg * 3e-3)
    for (fr in c(1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3)) {
      v <- frm_curve_deriv(fg, var = "x", order = o, newdata = nd, simultaneous = FALSE, eps = rg * fr)
      cat(sprintf("NCGPE design %s order %d step %.0e*range: value error %.2e (max |d| %.3g), se vs 3e-3: %.2e\n",
                  dn, o, fr, max(abs(v$.estimate - ex)), max(abs(ex)),
                  max(abs(v$.se / ref$.se - 1))))
    }
  }
}
