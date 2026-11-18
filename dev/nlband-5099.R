# Lane nlband, item 5099: the denominators the pole restart does not read,
# `EC50^h + x^h` (Hill, the limit battery's hill and fn_hillstep) and
# `1 + exp(...)` (logis, the fresh rows' logis3 and sig4). For each fit:
# nl_pole_walls() of every coefficient, and each denominator over a grid of
# each coefficient (the others at the fit) from -100 to 100 times its
# estimate's size: the smallest finite value, how many grid points give a
# non-finite one, and whether any row's denominator changes sign.
# NLBAND_ARM=after Rscript dev/nlband-5099.R
source("dev/nlband-common.R")
fits <- list()
set.seed(1)
d <- data.frame(x = runif(200, 0.1, 10))
d$y <- 10 * d$x^2 / (3^2 + d$x^2) + rnorm(200, 0, 0.5)
fits$hill <- frm(bf(y ~ Emax * x^h / (EC50^h + x^h), Emax + EC50 + h ~ 1,
                    nl = TRUE), data = d,
                 start = list(beta = c(10, 3, 2)))
set.seed(1)
d <- data.frame(x = c(runif(100, 0.5, 4.5), runif(100, 5.5, 10)))
d$y <- 10 * (d$x > 5) + rnorm(200, 0, 0.5)
fits$hillstep <- frm(bf(y ~ Emax * x^h / (EC50^h + x^h),
                        Emax + EC50 + h ~ 1, nl = TRUE), data = d,
                     start = list(beta = c(10, 5, 2)))
set.seed(1)
d <- data.frame(x = runif(200, 0, 10))
d$y <- 10 / (1 + exp((5 - d$x) / 1)) + rnorm(200, 0, 0.5)
fits$logis <- frm(bf(y ~ Asym / (1 + exp((xmid - x) / scal)),
                     Asym + xmid + scal ~ 1, nl = TRUE), data = d,
                  start = list(beta = c(10, 5, 1)))
fits$sig4 <- frm(bf(y ~ lo + (hi - lo) / (1 + exp((xmid - x) / scal)),
                    lo + hi + xmid + scal ~ 1, nl = TRUE), data = d,
                 start = list(beta = c(0, 10, 5, 1)))
for (nm in names(fits)) {
  f <- suppressWarnings(fits[[nm]])
  lp <- Filter(function(q) !is.null(q$nl_body), f$frame$linpreds)[[1L]]
  dens <- ns$nl_denoms(lp$nl_body)
  b <- f$estimates$beta
  cat("FIT", nm, "logLik", format(llv(f), digits = 10), "beta",
      paste(signif(b, 5), collapse = " "), "dens",
      paste(vapply(dens, deparse1, ""), collapse = " ; "), "\n")
  for (at in seq_along(b)) {
    wl <- ns$nl_pole_walls(f$frame, lp, f$estimates, "beta", at, dens)
    sz <- max(1, abs(b[[at]]))
    grid <- seq(-100, 100, length.out = 4001) * sz
    mn <- Inf
    nonfin <- 0L
    flips <- 0L
    for (v in grid) {
      t <- f$estimates
      t$beta[at] <- v
      vals <- suppressWarnings(ns$nl_body_at(f$frame, lp, t, exprs = dens))
      for (k in seq_along(vals)) {
        a <- vals[[k]]
        if (!all(is.finite(a))) nonfin <- nonfin + 1L
        a <- a[is.finite(a)]
        if (length(a)) {
          mn <- min(mn, a)
          if (any(a < 0)) flips <- flips + 1L
        }
      }
    }
    cat(sprintf("  coef %d (%s): walls %s | grid %g to %g: min finite %.4g, non-finite at %d, negative at %d\n",
                at, names(b)[at], if (is.null(wl)) "NULL" else "READ",
                min(grid), max(grid), mn, nonfin, flips))
  }
}
cat("DONE\n")
