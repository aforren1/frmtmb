# Every pair of a survey rerun against dev/rel079-equiv-out: verdict, level, agree and
# the numbers a verdict rests on. Prints each pair that moves, and the count.
#   Rscript dev/rel079-equiv-diff.R <old folder> <new folder>
a <- commandArgs(TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character", na.strings = "NA")
ids <- sort(sub("[.]tsv$", "", list.files(a[1], "^p[0-9]+[.]tsv$")))
idn <- sort(sub("[.]tsv$", "", list.files(a[2], "^p[0-9]+[.]tsv$")))
cat("pairs", length(ids), "and", length(idn), "; only in one:", setdiff(union(ids, idn), intersect(ids, idn)), "\n")
libs <- vapply(idn, function(id) {
  l <- grep("^lib:", readLines(file.path(a[2], "logs", paste0(id, ".log")), warn = FALSE), value = TRUE)
  if (length(l)) l[1L] else "none"
}, "")
print(table(sub(" [0-9.]+$", "", libs)))
cat("FAILED logs:", sum(vapply(idn, function(id) {
  any(grepl("FAILED", readLines(file.path(a[2], "logs", paste0(id, ".log")), warn = FALSE)))
}, NA)), "\n")
cols <- c("verdict", "level", "agree", "frm_ok", "brms_ok", "l1_ok", "l2_ok", "l4_ok", "proposal")
nums <- c("l1_const", "l1_grad", "l2_metric", "l4_z")
nv <- 0L
nm <- 0L
for (id in intersect(ids, idn)) {
  x <- rd(file.path(a[1], paste0(id, ".tsv")))
  y <- rd(file.path(a[2], paste0(id, ".tsv")))
  dc <- cols[vapply(cols, function(k) !identical(x[[k]], y[[k]]), NA)]
  dn <- nums[vapply(nums, function(k) {
    u <- suppressWarnings(as.numeric(strsplit(as.character(x[[k]]), ",")[[1L]]))
    v <- suppressWarnings(as.numeric(strsplit(as.character(y[[k]]), ",")[[1L]]))
    if (length(u) != length(v)) return(TRUE)
    if (!length(u)) return(FALSE)
    any(is.na(u) != is.na(v)) || isTRUE(any(abs(u - v) > 1e-6 * pmax(1, abs(u)), na.rm = TRUE))
  }, NA)]
  if (length(dc)) nv <- nv + 1L
  if (length(dn)) nm <- nm + 1L
  if (length(dc) || length(dn)) {
    cat(id, x$feature_a, "x", x$feature_b, ":", x$verdict, x$level, "->", y$verdict, y$level, "| moved:",
        paste(c(dc, dn), collapse = ","), "\n")
    for (k in dn) cat("   ", k, x[[k]], "->", y[[k]], "\n")
  }
}
cat("VERDICT COLUMNS MOVED on", nv, "pairs; numbers moved on", nm, "pairs\n")
