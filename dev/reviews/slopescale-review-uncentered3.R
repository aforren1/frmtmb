# Review of lane slopescale: which uncentered slope columns leave the
# lane's default short? x = c0 + s z for several (c0, s), y ~ 1 + (1 + x
# | g) and y ~ x + (1 + x | g), seeds 1..5. The model is invariant to an
# affine change of x, so the fit on the centered, standardized column is
# the reference. Prints the lane default's shortfall, whether it
# engaged, and autoscale = TRUE's shortfall. Lane library only.
# Usage: Rscript slopescale-review-uncentered3.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
mk <- function(seed, s, c0) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  z <- rnorm(300)
  eta <- 0.5 + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.5)[g] * z + 0.3 * z
  data.frame(g = g, x = c0 + s * z, y = eta + rnorm(300))
}
q <- function(expr) {
  nw <- 0L
  f <- withCallingHandlers(expr, warning = function(w) {
    nw <<- nw + 1L
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  attr(f, "nw") <- nw
  f
}
rows <- list()
for (cc in list(c(5, 1), c(5, 0.1), c(2000, 3), c(1, 0.05), c(0.5, 0.02)))
  for (fo in c("y ~ 1 + (1 + x | g)", "y ~ x + (1 + x | g)"))
    for (seed in 1:5) {
      d <- mk(seed, cc[2], cc[1])
      fml <- stats::as.formula(fo)
      dc <- d
      dc$x <- (d$x - mean(d$x)) / sd(d$x)
      ref <- as.numeric(logLik(q(frm(fml, data = dc))))
      fd <- q(frm(fml, data = d))
      ft <- q(frm(fml, data = d, control = frmtmb_control(autoscale = TRUE)))
      rows[[length(rows) + 1L]] <- data.frame(
        c0 = cc[1], s = cc[2], form = fo, seed = seed,
        engaged = isTRUE(fd$autoscaled),
        short_default = signif(ref - as.numeric(logLik(fd)), 3),
        warn_default = attr(fd, "nw"),
        short_TRUE = signif(ref - as.numeric(logLik(ft)), 3))
    }
r <- do.call(rbind, rows)
options(width = 200)
print(r, row.names = FALSE)
agg <- aggregate(cbind(short = r$short_default > 1e-3,
                       silent = r$short_default > 1e-3 & r$warn_default == 0,
                       shortT = r$short_TRUE > 1e-3) ~ c0 + s + form,
                 data = r, FUN = sum)
print(agg, row.names = FALSE)
