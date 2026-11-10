# Reviewer scratch (fix-round check of lane translator): the audit's
# coverage against the record's L1 pairs, and the blocks whose SD vector
# holds a small entry next to a large one (a half-collapsed block that a
# whole-block reversal would call seen). Reads files only.
ids <- scan("dev/translator-audit-ids.txt", what = "", quiet = TRUE)
r <- read.delim("dev/equiv-data/results.tsv", colClasses = "character")
l1 <- unique(r$id[grepl("L1", r$level)])
cat("audit ids", length(unique(ids)), "; record pairs with L1 in level",
    length(l1), "\n")
cat("L1 in record, not audited:", setdiff(l1, ids), "\n")
cat("audited, not L1 in record:", setdiff(ids, l1), "\n")
fs <- list.files("dev/translator-audit-out", "^p.*tsv$", full.names = TRUE)
a <- do.call(rbind, lapply(fs, function(p) {
  d <- read.delim(p, colClasses = "character")
  d[, c("id", "block", "max_delta", "tol", "blind", "sds")]
}))
a <- a[nzchar(a$block), ]
a$ratio <- as.numeric(a$max_delta) / as.numeric(a$tol)
num <- function(s) {
  v <- unlist(regmatches(s, gregexpr("[-0-9.e]+(?=(/| |$))", s,
                                     perl = TRUE)))
  suppressWarnings(as.numeric(v))
}
a$min_sd <- vapply(a$sds, function(s) suppressWarnings(min(num(s),
                                                            na.rm = TRUE)), 0)
op <- options(width = 250)
cat("\nblocks in pairs with some SD below 1e-2:\n")
print(a[is.finite(a$min_sd) & a$min_sd < 1e-2,
        c("id", "block", "ratio", "blind", "sds")], row.names = FALSE)
cat("\nlowest ratios max_delta / tol:\n")
print(head(a[order(a$ratio), c("id", "block", "ratio", "sds")], 12),
      row.names = FALSE)
options(op)
