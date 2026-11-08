# Item 1 summary: the default fit of each arm against the reference
# (dev/fitedge-gpbody-sweep.R). Short means more than 1e-3 nats below the
# reference; silent means conv 0 and no warning.
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
rd <- function(f) read.delim(file.path(o, f), stringsAsFactors = FALSE)
ref <- rd("gpbody-ref-lane.tsv")
k <- c("body", "curve", "n", "seed")
for (arm in commandArgs(TRUE)) {
  fit <- rd(sprintf("gpbody-fit-%s.tsv", arm))
  m <- merge(fit, ref[, c(k, "ll")], by = k, suffixes = c("", ".ref"))
  m$gap <- m$ll - pmax(m$ll.ref, m$ll)
  m$above <- m$ll - m$ll.ref
  m$short <- m$gap < -1e-3 | is.na(m$ll)
  m$silent <- m$short & m$conv %in% 0 & m$nwarn %in% 0
  cat(sprintf("SUM %s rows %d short %d silent %d worst %.3f above_ref_max %.2e errors %d time %.1f s\n",
              arm, nrow(m), sum(m$short), sum(m$silent),
              min(m$gap, na.rm = TRUE), max(m$above, na.rm = TRUE),
              sum(is.na(m$ll)), sum(m$sec)))
  for (b in unique(m$body)) {
    s <- m[m$body == b, ]
    cat(sprintf("SUM %s body %s rows %d short %d silent %d worst %.3f\n", arm,
                b, nrow(s), sum(s$short), sum(s$silent), min(s$gap, na.rm = TRUE)))
  }
  sh <- m[m$short, ]
  for (i in seq_len(nrow(sh))) {
    cat(sprintf("ROW %s %s %s n %d seed %d ll %.4f ref %.4f gap %.3f conv %s w %s %s\n",
                arm, sh$body[i], sh$curve[i], sh$n[i], sh$seed[i], sh$ll[i],
                sh$ll.ref[i], sh$gap[i], sh$conv[i], sh$nwarn[i], sh$wclass[i]))
  }
}
