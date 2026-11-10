# Compare two directories of per-pair tsvs (dev/equiv-measure.R output)
# column by column. `seconds` is wall time and differs on every run; the
# edited instrument's two new columns are listed, not compared.
#   Rscript dev/compat-noreg-compare.R <dir A> <dir B>
args <- commandArgs(trailingOnly = TRUE)
rd <- function(p) utils::read.delim(p, stringsAsFactors = FALSE,
                                    colClasses = "character")
fa <- sort(list.files(args[1], "[.]tsv$"))
fb <- sort(list.files(args[2], "[.]tsv$"))
ids <- intersect(fa, fb)
cat("pairs in A:", length(fa), "; in B:", length(fb), "; both:",
    length(ids), "\n")
new_cols <- character()
nd <- 0L
for (f in ids) {
  a <- rd(file.path(args[1], f))
  b <- rd(file.path(args[2], f))
  new_cols <- union(new_cols, setdiff(names(b), names(a)))
  miss <- setdiff(names(a), names(b))
  cols <- setdiff(intersect(names(a), names(b)), "seconds")
  diff <- cols[a[1, cols] != b[1, cols] |
                 xor(is.na(a[1, cols]), is.na(b[1, cols]))]
  diff <- diff[!is.na(diff)]
  st <- if ("stratum_status" %in% names(b))
    paste(b$stratum_status, "/", b$live_status) else ""
  cat(sprintf("%-6s %-26s %-16s %s\n", sub("[.]tsv$", "", f), a$verdict,
              st, if (length(diff) || length(miss)) paste(
                "DIFFERS:", paste(c(diff, if (length(miss))
                  paste0("-", miss)),
                                  collapse = " ")) else "identical"))
  for (k in diff) {
    cat("    ", k, ": ", substr(a[[k]], 1, 80), " | ", substr(b[[k]], 1, 80),
        "\n", sep = "")
  }
  if (length(diff) || length(miss)) nd <- nd + 1L
}
cat("columns only in B:", paste(new_cols, collapse = " "), "\n")
cat("pairs with any difference outside seconds:", nd, "of", length(ids),
    "\n")
