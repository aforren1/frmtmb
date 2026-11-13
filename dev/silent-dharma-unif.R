# Item 650: are dharma_residuals() uniform on a correctly specified
# cens() fit, and do they still reject a wrong one? 20 seeds each.
# Usage: Rscript dev/silent-dharma-unif.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
gauss_data <- function(seed, n = 300, skew = FALSE) {
  set.seed(seed)
  x <- rnorm(n)
  e <- if (skew) (rexp(n) - 1) else rnorm(n)
  ys <- 1 + 0.5 * x + e
  rc <- runif(n, 0.8, 2.5)
  lc <- runif(n, -1, 0.3)
  cen <- ifelse(ys > rc, 1, ifelse(ys < lc, -1, 0))
  y <- ifelse(cen == 1, rc, ifelse(cen == -1, lc, ys))
  hi <- y
  iv <- which(cen == 0)[seq(1, sum(cen == 0), by = 6)]
  cen[iv] <- 2
  y[iv] <- floor(ys[iv] * 2) / 2
  hi[iv] <- y[iv] + 0.5
  data.frame(x = x, y = y, cen = cen, hi = hi)
}
pois_data <- function(seed, n = 300) {
  set.seed(seed)
  x <- rnorm(n)
  ys <- rpois(n, exp(1 + 0.4 * x))
  k <- sample(1:4, n, TRUE)
  cen <- as.numeric(ys > k)
  data.frame(x = x, y = ifelse(cen == 1, k, ys), cen = cen)
}
ks_p <- function(r) suppressWarnings(stats::ks.test(r, "punif")$p.value)
out <- NULL
for (s in 1:20) {
  dg <- gauss_data(s)
  fg <- frm(bf(y | cens(cen, hi) ~ x), family = gaussian(), data = dg)
  rg <- dharma_residuals(fg, nsim = 250, seed = s)$scaledResiduals
  dw <- gauss_data(s, skew = TRUE)
  fw <- frm(bf(y | cens(cen, hi) ~ x), family = gaussian(), data = dw)
  rw <- dharma_residuals(fw, nsim = 250, seed = s)$scaledResiduals
  dp <- pois_data(s)
  fp <- frm(bf(y | cens(cen) ~ x), family = poisson(), data = dp)
  rp <- dharma_residuals(fp, nsim = 250, seed = s)$scaledResiduals
  out <- rbind(out, data.frame(
    seed = s, cens_share_g = mean(dg$cen != 0),
    p_gauss = ks_p(rg), p_wrong = ks_p(rw), p_pois = ks_p(rp),
    cens_share_p = mean(dp$cen != 0)))
}
print(format(out, digits = 3), row.names = FALSE)
cat(sprintf("SUMMARY %s: p<0.01 on gaussian %d/20, on poisson %d/20, on the wrong (skewed) model %d/20; median p gaussian %.3g poisson %.3g wrong %.3g\n",
            arm, sum(out$p_gauss < 0.01), sum(out$p_pois < 0.01),
            sum(out$p_wrong < 0.01), median(out$p_gauss),
            median(out$p_pois), median(out$p_wrong)))
