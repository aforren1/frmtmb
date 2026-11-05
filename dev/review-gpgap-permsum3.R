# Short check of fix round 1, item 1: Rhat and ESS over each set of one-chain runs.
#   Rscript dev/review-gpgap-permsum3.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (set in c("1?", "3?")) {
  for (arm in c("same", "perm")) {
    fs <- Sys.glob(sprintf("dev/review-gpgap-out/permdraw-%s-15-%s.rds", arm, set))
    cat(sprintf("set %s %s: %d chains\n", set, arm, length(fs)))
    if (length(fs) < 2) next
    ch <- lapply(fs, readRDS)
    for (v in names(ch[[1]])) {
      m <- do.call(cbind, lapply(ch, function(d) d[[v]][seq_len(min(vapply(ch, nrow, 0L)))]))
      q <- quantile(as.numeric(m), c(.05, .5, .95))
      cat(sprintf("PSUM %s %s %-11s q05 %.4g q50 %.4g q95 %.4g rhat %.3f bulk ess %.0f\n", set,
                  arm, v, q[1], q[2], q[3], posterior::rhat(m), posterior::ess_bulk(m)))
    }
  }
}
