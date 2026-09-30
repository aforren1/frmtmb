wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
l <- readRDS(file.path(wt, "dev/ordinal-rev-bitwise-lane.rds"))
for (m in c("cum_logit", "hurdle", "mv_ord_gauss", "categorical", "gamma")) {
  for (k in c("emm", "resid", "loglik_rows", "hyp", "predict", "sim"))
    if (is.character(l[[m]][[k]]) && length(l[[m]][[k]]) == 1)
      cat(m, k, ":", substr(l[[m]][[k]], 1, 160), "\n")
}
