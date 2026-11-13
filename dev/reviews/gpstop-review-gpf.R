# Review of lane gpstop: fitedge's 192 rows recounted from the lane's own
# TSVs (dev/gpstop-out/gpf-<set>-<arm>.tsv), independently of gpstop-gpsum.R.
setwd("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-out")
for (s in c("fresh", "lane96")) {
  r <- lapply(c("base", "lane", "ref"), function(a) {
    read.delim(sprintf("gpf-%s-%s.tsv", s, a))
  })
  k <- paste(r[[1]]$design, r[[1]]$curve, r[[1]]$n, r[[1]]$seed)
  stopifnot(identical(k, paste(r[[2]]$design, r[[2]]$curve, r[[2]]$n,
                               r[[2]]$seed)))
  best <- pmax(r[[1]]$ll, r[[2]]$ll, r[[3]]$ll, na.rm = TRUE)
  d <- r[[2]]$ll - r[[1]]$ll
  cat(s, "rows", length(d), "up>1e-6", sum(d > 1e-6), "down<-1e-6",
      sum(d < -1e-6), "max up", max(d), "min", min(d), "\n")
  cat(s, "short base", sum(best - r[[1]]$ll > 1e-3), "lane",
      sum(best - r[[2]]$ll > 1e-3), "\n")
  cat(s, "conv!=0 base", sum(r[[1]]$conv != 0), "lane", sum(r[[2]]$conv != 0),
      "| nwarn>0 base", sum(r[[1]]$nwarn > 0), "lane", sum(r[[2]]$nwarn > 0),
      "\n")
  i <- which(abs(d) > 1e-6 | r[[1]]$nwarn != r[[2]]$nwarn |
               r[[1]]$conv != r[[2]]$conv)
  if (length(i)) {
    print(cbind(r[[2]][i, 1:4], base = r[[1]]$ll[i], lane = r[[2]]$ll[i],
                ref = r[[3]]$ll[i], wb = r[[1]]$nwarn[i], wl = r[[2]]$nwarn[i],
                th_b = r[[1]]$theta[i], th_l = r[[2]]$theta[i]))
  }
  cat(s, "time base", sum(r[[1]]$sec), "lane", sum(r[[2]]$sec), "\n")
}
