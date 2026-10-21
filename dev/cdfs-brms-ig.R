# brms 2.23.0's inv_gaussian_lccdf() (inst/chunks/fun_inv_gaussian.stan)
# written out in R, against the survivor 1 - F and statmod's upper tail
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
brms_lccdf <- function(y, mu, shape) {
  log1p(-(pnorm(sqrt(shape) / sqrt(y) * (y / mu - 1)) -
            exp(2 * shape / mu) * pnorm(-sqrt(shape) / sqrt(y) * (y / mu + 1))))
}
brms_lcdf <- function(y, mu, shape) {
  log(pnorm(sqrt(shape) / sqrt(y) * (y / mu - 1)) +
        exp(2 * shape / mu) * pnorm(-sqrt(shape) / sqrt(y) * (y / mu + 1)))
}
for (y in c(0.5, 1, 2, 4)) {
  cat(sprintf("y %.1f mu 2 shape 4: brms lccdf %.10f  log(1 - exp(brms lcdf)) %.10f  statmod %.10f\n",
              y, brms_lccdf(y, 2, 4), log1p(-exp(brms_lcdf(y, 2, 4))),
              statmod::pinvgauss(y, 2, 4, lower.tail = FALSE, log.p = TRUE)))
}
