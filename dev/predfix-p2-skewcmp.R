# Per-fit comparison of the skew_normal set between the round-1 default
# (always seeded from the pre-fit) and the round-2 default (falls back
# when the pre-fit's coefficient block is flat or singular).
rd <- function(f) {
  x <- readLines(f)
  x <- grep("^(alpha0|alpha1|mu) ", x, value = TRUE)
  do.call(rbind, lapply(strsplit(trimws(x), " +"), function(v)
    data.frame(design = v[1], seed = as.integer(v[2]), s = as.numeric(v[3]),
               def = as.numeric(v[4]), false = as.numeric(v[5]))))
}
d <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
a <- rd(paste0(d, "p1-skewdown.txt"))
b <- rd(paste0(d, "p2-rerun-p1-skewdown.txt"))
m <- merge(a, b, by = c("design", "seed", "s"), suffixes = c(".r1", ".r2"))
m$lost <- m$def.r1 - m$def.r2
cat("fits", nrow(m), "; round-2 default below round-1 by more than 1e-6:",
    sum(m$lost > 1e-6), "; largest loss", signif(max(m$lost), 4),
    "; round-2 below FALSE by more than 1e-6:", sum(m$def.r2 - m$false.r2 < -1e-6), "\n")
print(m[m$lost > 1e-6, c("design", "seed", "s", "def.r1", "def.r2", "false.r2", "lost")],
      digits = 10)
