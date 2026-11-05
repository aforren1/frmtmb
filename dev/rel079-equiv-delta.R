# Every verdict change of the recorded 0.79.0 survey against the 0.78.0 baseline,
# and the corrected-weights line under the two frames.
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
a <- rd("C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out/results.tsv")
b <- rd("dev/equiv-data/results.tsv")
m <- merge(a[, c("id", "feature_a", "feature_b", "status", "verdict", "level", "agree")],
           b[, c("id", "verdict", "level", "agree", "reason")], by = "id", suffixes = c("_old", "_new"))
ch <- m[m$verdict_old != m$verdict_new | m$level_old %in% c(NA) != m$level_new %in% c(NA) |
          (!is.na(m$level_old) & !is.na(m$level_new) & m$level_old != m$level_new) |
          m$agree_old != m$agree_new, ]
cat("pairs", nrow(m), " changed", nrow(ch), "\n")
for (i in seq_len(nrow(ch))) with(ch[i, ], cat(sprintf(
  "%s %s x %s [%s]: %s %s (%s) -> %s %s (%s)\n", id, feature_a, feature_b, status, verdict_old, level_old,
  agree_old, verdict_new, level_new, agree_new)))
print(table(b$verdict))
fr <- utils::read.delim("dev/equiv-frame-out/compat.tsv", colClasses = "character")
r17 <- utils::read.delim("dev/rel079-log/frame-r17/compat.tsv", colClasses = "character")
k <- function(x) paste(x$feature_a, x$feature_b)
cat("core pairs r15", sum(fr$in_core_only == "TRUE"), " r17", sum(r17$in_core_only == "TRUE"), "\n")
print(table(fr$status[fr$in_core_only == "TRUE"])); print(table(r17$status[r17$in_core_only == "TRUE"]))
m2 <- merge(fr[fr$in_core_only == "TRUE", c("feature_a","feature_b","status")], r17[r17$in_core_only == "TRUE", c("feature_a","feature_b","status")], by = c("feature_a","feature_b"))
cat("core pairs whose status moved:", sum(m2$status.x != m2$status.y), "of", nrow(m2), "\n")
print(table(paste(m2$status.x, "->", m2$status.y)[m2$status.x != m2$status.y]))
