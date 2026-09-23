# Cross-check the gap range the wt-drmtmb lane recorded for the stated
# stream (6.406930 to 35.137329) against this lane's own run.
d <- read.csv(paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/",
                     "dev/skewinit-log/starts.csv"))
for (st in unique(d$stream)) {
  s <- d[d$stream == st, ]
  bad <- s$default < -1e-6
  cat(sprintf("%-10s stalls %2d of %d   gaps %.6f to %.6f\n", st, sum(bad),
              nrow(s), -max(s$default[bad]), -min(s$default[bad])))
}
