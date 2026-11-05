# Verdict changes of dev/refuse-equiv-out against the 0.78.0 baseline.
b <- read.delim("C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out/results.tsv")
m <- read.delim("dev/refuse-equiv-out/results.tsv")
x <- merge(b[, c("id", "feature_a", "feature_b", "verdict", "level")],
           m[, c("id", "verdict", "level", "l1_const", "l1_grad", "reason")],
           by = "id", suffixes = c("_base", "_lane"))
cat("pairs:", nrow(x), "\n")
print(table(base = x$verdict_base, lane = x$verdict_lane))
ch <- x[x$verdict_base != x$verdict_lane |
          !mapply(identical, as.character(x$level_base),
                  as.character(x$level_lane)), ]
cat("\nchanged verdict or level:\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf("%s %s x %s: %s (%s) -> %s (%s); l1_const %s grad %s\n",
              ch$id[i], ch$feature_a[i], ch$feature_b[i], ch$verdict_base[i],
              ch$level_base[i], ch$verdict_lane[i], ch$level_lane[i],
              format(ch$l1_const[i], digits = 3),
              format(ch$l1_grad[i], digits = 3)))
}
d <- ch[ch$verdict_lane %in% c("defect", "cannot measure") &
          ch$verdict_base != ch$verdict_lane, ]
if (nrow(d)) {
  cat("\nnew defect or cannot measure:\n")
  for (i in seq_len(nrow(d))) cat(d$id[i], substr(d$reason[i], 1, 300), "\n")
}
s <- readLines("dev/refuse-equiv-out/summary.md")
cat("\n", grep("headline", s, value = TRUE), "\n")
s0 <- readLines("C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out/summary.md")
cat("base:", grep("headline", s0, value = TRUE), "\n")
