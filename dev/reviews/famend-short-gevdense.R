# Re-check of lane famend after fix round 2: a dense xi scan of outside
# GEV rows, several mu and sigma, counting non-finite gradients of lpdf
# (and of lcdf and lccdf on a subset), and NaN Hessian entries on every
# 7th row.
#   Rscript dev/reviews/famend-short-gevdense.R <lib> <part 1|2|3>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
cat("lib", dirname(find.package("frmtmb")), "\n")
fam <- gen_extreme_value()
ev <- function(slot, y, xi, mu, ls, hess = FALSE) {
  p <- c(mu, ls, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(g = tp$jacobian(p), h = if (hess) tp$jacfun()$jacobian(p))
}
part <- as.integer(a[2])
xis <- switch(part,
  seq(0.40, 0.90, by = 0.005),
  c(-0.9, -0.7, -0.5, -0.3, -0.1, -0.03, -0.01, -1e-3, -1e-4, -1e-8,
    1e-8, 1e-4, 1e-3, 0.01, 0.03, 0.06, 0.1, 0.2, 0.3, 0.35, 0.95, 1, 1.5,
    2, 3, 5, 7, 10),
  c(0.41, 0.47, 0.53, 0.61, 0.66, 0.73, 0.79, 0.83, 0.86))
ms <- list(c(0, 0), c(3, 0.7), c(-2, -1.5), c(100, 2))
if (part == 3) ms <- list(c(0.3, -0.2), c(-50, 4), c(1e4, 0.1))
set.seed(part)
tot <- 0
bad <- 0
badh <- 0
nh <- 0
badc <- 0
nc <- 0
ex <- character()
for (xi in xis) {
  e <- if (xi > 0) -1 / xi else -1 / xi
  # outside offsets in scale units past the support's end: dense over the
  # band of h in (-700, -650) and spread to far away
  off <- c(seq(1e-6, 60, length.out = 120), 10^stats::runif(40, -10, 3),
           if (xi > 0) (seq(19, 21, by = 0.01) - 1 / xi)[
             seq(19, 21, by = 0.01) - 1 / xi > 0],
           1e4, 1e8, 1e12, 1e16)
  sgn <- if (xi > 0) -1 else 1
  for (m in ms) {
    k <- 0
    for (o in off) {
      x <- e + sgn * o
      if (!(xi * x < -1)) next
      y <- m[1] + exp(m[2]) * x
      k <- k + 1
      hh <- (k %% 7 == 0)
      r <- ev("lpdf", y, xi, m[1], m[2], hess = hh)
      tot <- tot + 1
      if (any(!is.finite(r$g))) {
        bad <- bad + 1
        if (length(ex) < 10) {
          ex <- c(ex, sprintf("xi %g mu %g ls %g x %.6g", xi, m[1], m[2], x))
        }
      }
      if (hh) {
        nh <- nh + 1
        badh <- badh + any(is.nan(r$h))
      }
      if (k %% 11 == 0) {
        for (s in c("lcdf", "lccdf")) {
          nc <- nc + 1
          badc <- badc + any(!is.finite(ev(s, y, xi, m[1], m[2])$g))
        }
      }
    }
  }
}
cat(sprintf(paste0("PART %d: %d xi x %d (mu, sigma); outside rows %d, ",
                   "lpdf non-finite gradient %d; NaN Hessian %d of %d; ",
                   "lcdf/lccdf non-finite gradient %d of %d\n"),
            part, length(xis), length(ms), tot, bad, badh, nh, badc, nc))
if (length(ex)) cat(ex, sep = "\n")
