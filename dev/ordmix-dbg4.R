.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d
}
d <- omx_data(20261011, n = 300)
print(table(d$y, d$g))
ob <- frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()), data = d, dry_run = "objective")
obj <- ob$obj
gr <- function(q) { g <- obj$gr(q); if (any(!is.finite(g))) { cat("NaN grad at\n"); print(q); print(g); print(obj$fn(q)); stop("nan") }; g }
o <- try(nlminb(obj$par, obj$fn, gr))
