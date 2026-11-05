# Rhat and bulk ESS across dev/ncgp-permsamp.R's chains, one row per arm and size.
#   Rscript dev/ncgp-permsamp-sum.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
fs <- list.files("dev/ncgp-out/permsamp", pattern = "[.]rds$", full.names = TRUE)
key <- sub("-[0-9]+[.]rds$", "", basename(fs))
for (k in unique(key)) {
  ch <- lapply(fs[key == k], readRDS)
  nd <- min(vapply(ch, nrow, 0L))
  for (v in names(ch[[1L]])) {
    m <- vapply(ch, function(d) d[[v]][seq_len(nd)], numeric(nd))
    cat(sprintf("NCGPPSS %-12s chains %d %-11s q05 %.3e q50 %.3e q95 %.3e rhat %.3f ess_bulk %.0f\n",
                k, length(ch), v, stats::quantile(m, 0.05), stats::median(m),
                stats::quantile(m, 0.95), posterior::rhat(m), posterior::ess_bulk(m)))
  }
}
