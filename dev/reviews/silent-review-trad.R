# dharma_residuals(method = "traditional" / rotation) on a censored fit
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5)
x <- rnorm(300); ys <- 1 + 0.5 * x + rnorm(300)
ct <- runif(300, 0.8, 2.5)
cen <- as.numeric(ys > ct)
d <- data.frame(x = x, y = ifelse(cen == 1, ct, ys), cen = cen)
f <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = d)
for (m in c("PIT", "traditional")) {
  r <- suppressMessages(dharma_residuals(f, nsim = 250, seed = 2,
                                         method = m))$scaledResiduals
  cat(sprintf("%s: right-censored rows %d, residual == 1 on %d; KS p %.3g\n",
              m, sum(cen == 1), sum(r[cen == 1] == 1),
              suppressWarnings(ks.test(r, "punif")$p.value)))
}
r <- suppressMessages(dharma_residuals(f, nsim = 250, seed = 2,
                                       rotation = "estimated"))
cat(sprintf("rotation: KS p %.3g\n",
            suppressWarnings(ks.test(r$scaledResiduals, "punif")$p.value)))
