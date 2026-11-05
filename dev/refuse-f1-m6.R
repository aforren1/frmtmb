# m6: a von_mises value a density reads, near the window's end: frm() holds it
# on the line now. Laplace against integrate() per value, on the review's
# shapes (dev/review-refuse-misc.R, seed 10201): mean 3.14 kappa 20, mean 3.05
# kappa 2, mean 0 kappa 2. Mean direction against the observed rows' fit.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(10201)
n <- 200
for (cfg in list(c(3.14, 20), c(3.05, 2), c(0, 2))) {
  a <- cfg[1]
  kap <- cfg[2]
  ang <- a + rnorm(n, 0, 1 / sqrt(kap))
  y <- ((ang + pi) %% (2 * pi)) - pi
  d <- data.frame(y = y)
  # the reader sees the value's cosine distance, so that it says nothing of
  # the period: linear in the value on its own branch
  d$z <- 0.2 * y + rnorm(n, 0, 0.5)
  d$y[seq(5, 195, 10)] <- NA
  fit <- suppressWarnings(frm(mvbf(bf(z ~ mi(y)),
                                   bf(y | mi() ~ 1, family = von_mises()),
                                   rescor = FALSE), data = d))
  env <- fit$obj$env
  pv <- env$last.par.best
  mi_i <- which(names(pv) == "miss")
  base <- env$f(pv)
  err <- numeric(length(mi_i))
  for (k in seq_along(mi_i)) {
    u0 <- pv[mi_i[k]]
    h <- function(u) {
      vapply(u, function(ui) exp(base - env$f(replace(pv, mi_i[k], ui))), 0)
    }
    d2 <- (env$f(replace(pv, mi_i[k], u0 + 1e-3)) +
             env$f(replace(pv, mi_i[k], u0 - 1e-3)) - 2 * base) / 1e-6
    w <- 1 / sqrt(max(d2, 1e-6))
    # one period, as the density is periodic in the value alone
    ex <- integrate(h, u0 - min(40 * w, pi), u0 + min(40 * w, pi),
                    subdivisions = 2000L, rel.tol = 1e-9,
                    stop.on.error = FALSE)$value
    err[k] <- log(sqrt(2 * pi) * w) - log(ex)
  }
  cc <- suppressWarnings(frm(bf(y ~ 1), family = von_mises(),
                             data = d[!is.na(d$y), ]))
  mj <- frmtmb:::fixef_by_dpar(fit)$y_mu[[1L]]
  mc <- frmtmb:::fixef_by_dpar(cc)$mu[[1L]]
  cat(sprintf(paste0("mean %.2f kappa %g: chart %s, Laplace - exact per value",
                     " %+.3f [%+.3f, %+.3f]; mean eta read %.3f complete %.3f;",
                     " conv %d\n"),
              a, kap, fit$frame$mi_map$y$chart, mean(err), min(err), max(err),
              mj, mc, fit$opt$convergence))
}
