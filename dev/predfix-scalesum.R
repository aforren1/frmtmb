# Summary numbers for item 6, generated from the saved calibration runs.
d <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
A <- readRDS(paste0(d, "scalecal-A.rds"))
cat("A: fits", nrow(A), "; default short > 1e-6 at scale 1e-3:",
    sum(A$def_bad[A$scale == 1e-3]), "of", sum(A$scale == 1e-3),
    ", max", signif(max(A$short_default[A$scale == 1e-3]), 3),
    "; at 1e-2:", sum(A$def_bad[A$scale == 1e-2]), "\n")
B <- readRDS(paste0(d, "scalecal-B.rds"))
for (m in c("ML", "REML", "profile")) {
  k <- B$mode == m
  cat("B", m, ": default short > 1e-6 in", sum(B$short_default[k] > 1e-6),
      "of", sum(k), ", max", signif(max(B$short_default[k]), 4),
      "; default equals autoscale within 1e-9 in",
      sum(abs(B$short_default[k] - B$short_auto[k]) < 1e-9), "\n")
}
r <- readRDS(paste0(d, "scalescan.rds"))
cat("scan: fits", nrow(r), "; scale > 1, default short:",
    sum(r$scale > 1 & (is.na(r$def_short) | r$def_short > 1e-6)), "of",
    sum(r$scale > 1), "; scale 1e-5 to 1, default short:",
    sum(r$scale >= 1e-5 & r$scale <= 1 & r$def_short > 1e-6), "of",
    sum(r$scale >= 1e-5 & r$scale <= 1),
    "; max default short:", signif(max(r$def_short, na.rm = TRUE), 5),
    "; autoscale short anywhere:", sum(r$auto_short > 1e-6, na.rm = TRUE), "\n")
