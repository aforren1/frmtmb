# Item 9: quantiles of each arm's chains, n 20, true sigma 1e-2, seed 6620.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
.libPaths(REV_USER)
for (arm in c("lane", "gpgap", "brms")) {
  fl <- list.files(REV_OUT, sprintf("^sampcmp-%s-[0-9]+[.]rds$", arm), full.names = TRUE)
  if (!length(fl)) next
  ms <- lapply(fl, readRDS)
  a <- array(NA_real_, c(nrow(ms[[1]]), length(ms), ncol(ms[[1]])), dimnames = list(NULL, NULL, colnames(ms[[1]])))
  for (i in seq_along(ms)) a[, i, ] <- ms[[i]][, colnames(ms[[1]])]
  for (v in dimnames(a)[[3]]) {
    x <- a[, , v]
    cat(sprintf("RSQ %-5s chains %d %-11s q05 %.4g q50 %.4g q95 %.4g rhat %.3f ess %.0f\n", arm, ncol(x), v, quantile(x, .05), median(x), quantile(x, .95),
                posterior::rhat(x), posterior::ess_bulk(x)))
  }
}
