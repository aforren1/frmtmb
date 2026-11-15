# Short check of lane famend: which entries are non-finite on inside GEV
# rows at xi near 0 (the same on 0.82.0, famend-short-gevclamp.R).
#   Rscript dev/reviews/famend-short-xi0.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
fam <- gen_extreme_value()
ev <- function(slot, y, xi) {
  p <- c(0.3, -0.2, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam[[slot]](y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                        .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  list(v = tp(p), g = tp$jacobian(p), h = tp$jacfun()$jacobian(p))
}
tab <- list()
for (xi in c(-1e-4, -1e-8, -1e-12, 0, 1e-12, 1e-8, 1e-4)) {
  lo <- max(if (xi > 0) -1 / xi else -60, -1e9)
  hi <- min(if (xi < 0) -1 / xi else 60, 1e9)
  xs <- c(seq(lo, hi, length.out = 201)[-c(1, 201)], -5, -1, -1e-3, 0, 1e-3,
          1, 5, 20, -20)
  if (xi > 0) xs <- c(xs, lo + c(1e-6, 1e-3, 1) / xi)
  if (xi < 0) xs <- c(xs, hi - c(1e-6, 1e-3, 1) / abs(xi))
  xs <- xs[xi * xs > -1]
  for (slot in c("lpdf", "lcdf", "lccdf")) {
    for (x in xs) {
      r <- ev(slot, 0.3 + exp(-0.2) * x, xi)
      key <- sprintf("%s v %s g %s h %s", slot,
                     paste(unique(r$v[!is.finite(r$v)]), collapse = ","),
                     paste(unique(r$g[!is.finite(r$g)]), collapse = ","),
                     paste(unique(r$h[!is.finite(r$h)]), collapse = ","))
      if (!all(is.finite(unlist(r)))) {
        tab[[key]] <- c(tab[[key]], x)
      }
    }
  }
}
for (k in names(tab)) {
  cat(sprintf("%-50s rows %3d, x from %.4g to %.4g\n", k, length(tab[[k]]),
              min(tab[[k]]), max(tab[[k]])))
}
