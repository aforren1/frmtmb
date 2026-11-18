# Lane squant: frm_bootstrap() (parametric: it draws from the fitted
# family) against the Hessian and sandwich standard errors on the homo
# design of dev/squant-se.R (n 400, quantile 0.9, gaussian errors), whose
# Monte Carlo sd of the slope is in dev/squant-out/se-summary.txt.
# Seeds 20261117 to 20261121, nsim 200.
# Run: Rscript dev/squant-boot.R > dev/squant-out/boot.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
for (s in 1:3) {
  set.seed(20261116 + s)
  n <- 400
  dd <- data.frame(x = runif(n), row = seq_len(n))
  dd$y <- 1 + 2 * dd$x + rnorm(n)
  fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(), data = dd)
  b <- frm_bootstrap(fit, nsim = 200, seed = s)

  m <- if (is.list(b) && !is.null(b$t)) b$t else do.call(rbind, b)
  bs <- apply(m, 2, sd)
  cat(sprintf("seed %d: slope SE hessian %.4f sandwich %.4f bootstrap %.4f\n",
              s, sqrt(vcov(fit)[2, 2]),
              sqrt(vcov_cluster(fit, ~ row, "CR0")[2, 2]),
              bs[grep("x", names(bs))[1]]))
}
