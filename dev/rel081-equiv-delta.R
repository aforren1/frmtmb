# Every verdict change of the recorded 0.81.0 survey against the 0.80.0
# record (dev/equiv-data/results.tsv at e0cce906), and the frames compared.
#   Rscript dev/rel081-equiv-delta.R > dev/rel081-log/equiv-delta.txt
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
old <- tempfile(fileext = ".tsv")
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
writeLines(system2("git", c("-c", "core.autocrlf=false", "show",
                            "e0cce906:dev/equiv-data/results.tsv"),
                   stdout = TRUE), old)
a <- rd(old)
b <- rd("dev/equiv-data/results.tsv")
m <- merge(a[, c("id", "feature_a", "feature_b", "status", "verdict", "level", "agree")],
           b[, c("id", "verdict", "level", "agree", "reason")], by = "id",
           suffixes = c("_old", "_new"))
ch <- m[m$verdict_old != m$verdict_new | is.na(m$level_old) != is.na(m$level_new) |
          (!is.na(m$level_old) & !is.na(m$level_new) & m$level_old != m$level_new) |
          m$agree_old != m$agree_new, ]
cat("pairs", nrow(m), " changed", nrow(ch), "\n")
for (i in seq_len(nrow(ch))) with(ch[i, ], cat(sprintf(
  "%s %s x %s [%s]: %s %s (%s) -> %s %s (%s) | %s\n", id, feature_a, feature_b, status,
  verdict_old, level_old, agree_old, verdict_new, level_new, agree_new,
  substr(reason, 1, 160))))
cat("\nverdicts, 0.80.0 record:\n"); print(table(a$verdict))
cat("verdicts, 0.81.0:\n"); print(table(b$verdict))
for (nm in c("l1_const", "l1_grad", "l2_metric", "l4_z")) {
  if (!nm %in% names(a) || !nm %in% names(b)) next
  x <- merge(a[, c("id", nm)], b[, c("id", nm)], by = "id")
  u <- suppressWarnings(as.numeric(x[[2]])); v <- suppressWarnings(as.numeric(x[[3]]))
  mv <- (is.na(u) != is.na(v)) | (!is.na(u) & !is.na(v) &
                                   abs(u - v) > 1e-6 * pmax(1, abs(u)))
  # two infinite values of one sign agree
  mv[is.na(mv)] <- !(u[is.na(mv)] == v[is.na(mv)]) %in% TRUE
  cat(sprintf("%s moved on %d pairs: %s\n", nm, sum(mv),
              paste(x$id[mv], collapse = " ")))
}
fr <- utils::read.delim("dev/equiv-frame-out/compat.tsv", colClasses = "character")
r20 <- utils::read.delim("dev/rel081-log/frame-r20/compat.tsv", colClasses = "character")
cat("\ncore pairs r17 frame", sum(fr$in_core_only == "TRUE"), " r20 frame",
    sum(r20$in_core_only == "TRUE"), "\n")
m2 <- merge(fr[fr$in_core_only == "TRUE", c("feature_a", "feature_b", "status")],
            r20[r20$in_core_only == "TRUE", c("feature_a", "feature_b", "status")],
            by = c("feature_a", "feature_b"))
cat("core pairs whose status moved:", sum(m2$status.x != m2$status.y), "of", nrow(m2), "\n")
print(table(paste(m2$status.x, "->", m2$status.y)[m2$status.x != m2$status.y]))
# the lane each changed pair belongs to
tr_ids <- c("p021", "p042", "p044", "p047", "p048", "p069", "p131", "p136",
            "p143", "p022", "p023", "p094", "p109")
co_ids <- c("p010", "p011", "p014", "p059", "p146", "p154", "p158", "p166",
            "p167", "p169", "p171")
lane <- ifelse(ch$id %in% tr_ids, "translator",
               ifelse(ch$id %in% co_ids, "compat", "NONE"))
cat("\nchanged pairs by lane:\n"); print(table(lane))
cat("changed pairs of no lane:", paste(ch$id[lane == "NONE"], collapse = " "), "\n")
# adjudicated columns, where results.tsv carries them
adj <- intersect(c("verdict_adj", "adjudicated", "adj_verdict"), names(b))
cat("adjudication columns in the new record:", paste(adj, collapse = " "), "\n")
