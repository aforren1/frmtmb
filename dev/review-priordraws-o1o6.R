# Reviewer: O1 b_Intercept 95% quantile of each side against the exact
# prior; O6 prior predictive with brms's non-finite draws dropped.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
set.seed(1); x <- rnorm(60)
sdI <- sqrt(4 + mean(x)^2); q95 <- 1 + qnorm(0.95) * sdI
for (s in c("brms", "frm")) {
  v <- readRDS(file.path(LOG, sprintf("only-%s-O1.rds", s)))$draws$b_Intercept
  m <- matrix(v, ncol = 4)
  cat(sprintf("O1 %s q95 %.3f exact %.3f mcse %.3f z %.2f\n", s, quantile(v, 0.95), q95,
              posterior::mcse_quantile(m, 0.95), (quantile(v, 0.95) - q95) / posterior::mcse_quantile(m, 0.95)))
}
b <- readRDS(file.path(LOG, "only-brms-O6.rds"))$yrep
f <- readRDS(file.path(LOG, "only-frm-O6.rds"))$yrep
for (j in 1:3) {
  bj <- b[is.finite(b[, j]), j]
  cat(sprintf("O6 yrep[,%d] brms finite %d min %.3f | frm min %.3f | KS p %.3f\n", j,
              length(bj), min(bj), min(f[, j]), suppressWarnings(ks.test(bj, f[, j])$p.value)))
}
