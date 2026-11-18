# Lane squant: frm_bootstrap() on asym_laplace() (the base build,
# rellib-r22), the design of dev/squant-boot.R. The Monte Carlo sd of
# rq()'s slope there is 0.2977 (dev/squant-se.R, 300 reps).
# Run: Rscript dev/squant-boot-ald.R > dev/squant-out/boot-ald.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
for (s in 1:3) {
  set.seed(20261116 + s)
  n <- 400
  dd <- data.frame(x = runif(n))
  dd$y <- 1 + 2 * dd$x + rnorm(n)
  fit <- suppressWarnings(frm(bf(y ~ x, quantile = 0.9),
                              family = asym_laplace(), data = dd))
  b <- suppressWarnings(frm_bootstrap(fit, nsim = 200, seed = s))
  r <- quantreg::rq(y ~ x, tau = 0.9, data = dd)
  cat(sprintf(paste0("seed %d: slope bootstrap SE %.4f (converged %d of",
                     " 200) | rq nid SE %.4f\n"), s, sd(b$t[, 2]),
              sum(b$converged), summary(r, se = "nid")$coefficients[2, 2]))
}
