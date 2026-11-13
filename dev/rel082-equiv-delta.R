# Every verdict change of the recorded 0.82.0 survey against the 0.81.0
# record (dev/equiv-data/results.tsv at e8abed0a), and the frames compared.
#   Rscript dev/rel082-equiv-delta.R > dev/rel082-log/equiv-delta.txt
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
old <- tempfile(fileext = ".tsv")
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
writeLines(system2("git", c("-c", "core.autocrlf=false", "show",
                            "e8abed0a:dev/equiv-data/results.tsv"),
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
cat("\nverdicts, 0.81.0 record:\n"); print(table(a$verdict))
cat("verdicts, 0.82.0:\n"); print(table(b$verdict))
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
# no lane but those three remeasured pairs this round (names: p139 and
# p129; ordstart: p131); arspat states no pair status changed, and its
# five frm_compat() moves (dev/rel082-log/frame-diff-r20-r21.txt) name no
# surveyed pair
for (id in c("p139", "p129", "p131")) {
  cat(id, "in the changed set:", id %in% ch$id, "\n")
}
cat("changed pairs:", paste(ch$id, collapse = " "), "\n")
