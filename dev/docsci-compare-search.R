# Why the committed core search.json is 2.2 MB and the new build's is
# 1.2 MB. Compares the two indexes entry by entry rather than guessing
# from the byte count.
one <- function(x) {
  if (is.null(x) || !length(x)) "" else as.character(x)[1]
}
load1 <- function(p) jsonlite::fromJSON(p, simplifyVector = FALSE)
a <- load1("docs/search.json")
b <- load1("dev/docsci-site/search.json")
cat("committed entries: ", length(a), "\n", sep = "")
cat("new build entries: ", length(b), "\n", sep = "")
cat("committed keys: ", paste(unique(unlist(lapply(a, names))),
                              collapse = " "), "\n", sep = "")
pa <- vapply(a, function(x) one(x$path), character(1))
pb <- vapply(b, function(x) one(x$path), character(1))
only_a <- setdiff(pa, pb)
cat("paths only in committed: ", length(only_a), "\n", sep = "")
cat("paths only in new:       ", length(setdiff(pb, pa)), "\n", sep = "")
cat("\ncount of committed-only paths by directory:\n")
print(sort(table(dirname(only_a)), decreasing = TRUE)[1:10])
cat("\nfirst 15 committed-only paths:\n")
cat(paste0("  ", head(only_a, 15)), sep = "\n")
cat("\nall paths only in the new index:\n")
if (length(setdiff(pb, pa))) {
  cat(paste0("  ", setdiff(pb, pa)), sep = "\n")
} else {
  cat("  (none)\n")
}
nch <- function(x) sum(nchar(unlist(lapply(x, function(e) unlist(e)))))
cat("\ncommitted characters: ", nch(a), "\n", sep = "")
cat("new characters:       ", nch(b), "\n", sep = "")
cat("characters in committed-only entries: ",
    nch(a[pa %in% only_a]), "\n", sep = "")
