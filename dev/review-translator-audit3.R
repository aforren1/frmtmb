# Reviewer scratch (second fix-round check of lane translator): block
# shapes in the final per-row audit, the weakest parts, and whether any
# audited block has a prior that reversal does not keep (rcar, zcar, Ymi).
# Reads files only.
fs <- list.files("dev/translator-audit3-out", "^p.*tsv$", full.names = TRUE)
a <- do.call(rbind, lapply(fs, function(p) {
  d <- read.delim(p, colClasses = "character")
  d[, intersect(c("id", "pair", "block", "max_delta", "tol", "blind", "sds",
                  "points"), names(d))]
}))
a <- a[nzchar(a$block), ]
a$ratio <- as.numeric(a$max_delta) / as.numeric(a$tol)
a$kind <- sub("_.*$", "", sub("\\[.*$", "", a$block))
cat("parts:", nrow(a), " pairs:", length(unique(a$id)), "\n")
print(table(a$kind))
cat("\nrow parts:", sum(grepl("row", a$block)), " by parts:",
    sum(grepl(", by ", a$block)), "\n")
cat("\nnon-iid-prior blocks (rcar, zcar, Ymi, Yl):\n")
print(a[grepl("^(rcar|zcar|Ymi|Yl)", a$block),
        c("id", "pair", "block", "ratio")], row.names = FALSE)
op <- options(width = 250)
cat("\nweakest 15 parts:\n")
print(head(a[order(a$ratio), c("id", "pair", "block", "ratio", "points",
                               "sds")], 15), row.names = FALSE)
cat("\nmatrix blocks with a single row tested (no [row] entries):\n")
mat <- unique(sub("\\[.*$", "", a$block[grepl("\\[row", a$block)]))
cat(length(mat), "blocks have row parts\n")
options(op)
