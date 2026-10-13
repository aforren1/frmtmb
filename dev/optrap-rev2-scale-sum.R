# Reviewer of lane optrap, re-check: dev/optrap-rev2-scale.R's lane runs
# on build 26 against round 0's lane and rellib-r7 runs
# (dev/optrap-rev-log/scale-*-{lane,r7}.tsv). Reference: the best
# unscaled log-likelihood any of the three arms or the external fit
# reached. Per scale: short / nonzero code / convergence-warned.
cw <- function(w) grepl("maximum absolute gradient|did not report", w)
for (cs in c("student_sigx", "exgauss", "mixture", "skewn", "huber",
             "gauss_sigx")) {
  rd <- function(f) read.delim(f, quote = "")
  A <- list(r7 = rd(sprintf("dev/optrap-rev-log/scale-%s-r7.tsv", cs)),
            round0 = rd(sprintf("dev/optrap-rev-log/scale-%s-lane.tsv", cs)),
            build26 = rd(sprintf("dev/optrap-rev2-log/scale-%s-lane.tsv", cs)))
  A <- lapply(A, function(x) x[x$seed %in% 1:20, ])
  all <- do.call(rbind, A)
  ref <- tapply(pmax(all$llu, all$ext, na.rm = TRUE), all$seed, max,
                na.rm = TRUE)
  for (nm in names(A)) {
    x <- A[[nm]]
    gap <- x$llu - ref[as.character(x$seed)]
    out <- sapply(sort(unique(x$scale)), function(sc) {
      i <- x$scale == sc
      sprintf("%d/%d/%d", sum(gap[i] < -1e-4, na.rm = TRUE),
              sum(x$code[i] != 0, na.rm = TRUE), sum(cw(x$warn[i])))
    })
    cat(sprintf("%-12s %-8s %s\n", cs, nm, paste(out, collapse = "  ")))
  }
}
