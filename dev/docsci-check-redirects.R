# The 105 core reference pages that exist in the committed docs/ and not
# in the new build. Are they orphans or deliberate redirects, and are
# they in the core search index?
a <- jsonlite::fromJSON("docs/search.json", simplifyVector = FALSE)
one <- function(x) if (is.null(x) || !length(x)) "" else as.character(x)[1]
pa <- vapply(a, function(x) one(x$path), character(1))

docs <- sort(basename(Sys.glob("docs/reference/*.html")))
new <- sort(basename(Sys.glob("dev/docsci-site/reference/*.html")))
only <- setdiff(docs, new)
cat("core reference pages: docs ", length(docs), ", new ", length(new),
    ", only in docs ", length(only), "\n", sep = "")

is_stub <- function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  grepl("http-equiv=\"refresh\"", txt, fixed = TRUE)
}
stub <- vapply(file.path("docs/reference", only), is_stub, logical(1))
cat("of those, meta-refresh redirect stubs: ", sum(stub), "\n", sep = "")
cat("not stubs: ", paste(only[!stub], collapse = " "), "\n", sep = "")

# Is any of them in the CORE search index? The path there is an absolute
# URL, so the test anchors on the core prefix and not on a substring that
# a subsite path would also match.
core <- "https://aforren1.github.io/frmtmb/reference/"
hit <- only[paste0(core, only) %in% pa]
cat("stubs present in the core search index: ", length(hit), "\n",
    sep = "")
if (length(hit)) cat("  ", paste(head(hit, 10), collapse = " "), "\n",
                     sep = "")

# Where do the stubs point?
tgt <- vapply(file.path("docs/reference", only[stub]), function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(txt, regexpr("URL=[^\"]+", txt))
  if (length(m)) sub("^URL=", "", m) else NA_character_
}, character(1))
cat("\nredirect targets by subsite:\n")
print(sort(table(sub("^(https://[^/]+/frmtmb/[^/]+)/.*$", "\\1", tgt)),
           decreasing = TRUE))
cat("\ntargets that do not exist in the new build:\n")
rel <- sub("^https://aforren1[.]github[.]io/frmtmb/", "", tgt)
bad <- rel[!file.exists(file.path("dev/docsci-site", rel))]
cat(if (length(bad)) paste0("  ", bad, collapse = "\n") else "  (none)",
    "\n")
