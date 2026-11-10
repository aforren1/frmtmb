# Every pair whose status or note differs between two frm_compat() dumps
# (dev/equiv-dump-compat.R output), core pairs and all-loaded pairs.
#   Rscript dev/compat-frame-diff.R <before/compat.tsv> <after/compat.tsv>
args <- commandArgs(trailingOnly = TRUE)
rd <- function(f) {
  x <- utils::read.delim(f, stringsAsFactors = FALSE)
  a <- pmin(x$feature_a, x$feature_b)
  b <- pmax(x$feature_a, x$feature_b)
  x$key <- paste(a, b, sep = " x ")
  x
}
b <- rd(args[1])
a <- rd(args[2])
cat("before:", nrow(b), "pairs,", sum(b$in_core_only), "core\n")
cat("after: ", nrow(a), "pairs,", sum(a$in_core_only), "core\n")
cat("pairs only before:", length(setdiff(b$key, a$key)),
    "; only after:", length(setdiff(a$key, b$key)), "\n")
m <- merge(b[c("key", "status", "note", "in_core_only")],
           a[c("key", "status", "note", "in_core_only")],
           by = "key", suffixes = c("_before", "_after"))
st <- m[m$status_before != m$status_after, ]
nt <- m[m$status_before == m$status_after & m$note_before != m$note_after, ]
cat("\nstatus moves:", nrow(st), "(core", sum(st$in_core_only_after), ")\n")
print(table(paste(st$status_before, "->", st$status_after)))
for (i in order(st$key)) {
  cat(sprintf("  %-34s %-9s -> %-9s core=%s\n", st$key[i], st$status_before[i],
              st$status_after[i], st$in_core_only_after[i]))
}
cat("\nnote-only changes:", nrow(nt), "(core", sum(nt$in_core_only_after),
    ")\n")
for (i in order(nt$key)) {
  cat(sprintf("  %-34s %s\n", nt$key[i], nt$status_after[i]))
}
