# Short check of fix round 1, item 1: Rhat and ESS over the one-chain runs of
# dev/review-gpgap-perm.R, the patched factor against the control.
#   Rscript dev/review-gpgap-permsum.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (arm in c("same", "perm")) {
  fs <- Sys.glob(sprintf("dev/review-gpgap-out/permdraw-%s-15-1?.rds", arm))
  if (length(fs) < 2) {
    cat(arm, ": chains done", length(fs), "\n")
    next
  }
  ch <- lapply(fs, readRDS)
  cat(sprintf("%s: %d chains\n", arm, length(ch)))
  for (v in names(ch[[1]])) {
    m <- do.call(cbind, lapply(ch, function(d) d[[v]][seq_len(min(vapply(ch, nrow, 0L)))]))
    q <- quantile(as.numeric(m), c(.05, .5, .95))
    cat(sprintf("PSUM %s %-11s q05 %.4g q50 %.4g q95 %.4g rhat %.3f bulk ess %.0f tail ess %.0f\n",
                arm, v, q[1], q[2], q[3], posterior::rhat(m), posterior::ess_bulk(m),
                posterior::ess_tail(m)))
  }
}
