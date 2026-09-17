# Reviewer, lane wt-conditions: every stopifnot() and match.arg() call in
# the R/ trees of all 8 packages, with its enclosing top-level function,
# whether that function is exported, and the argument it checks.
#   Rscript dev/conditions-rev-unclassed-sites.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
trees <- c(frmtmb = root, setNames(
  file.path(root, "extensions", list.files(file.path(root, "extensions"))),
  list.files(file.path(root, "extensions"))))
rows <- list()
for (p in names(trees)) {
  exports <- getNamespaceExports(p)
  s3 <- unclass(asNamespace(p)[[".__S3MethodsTable__."]])
  for (f in list.files(file.path(trees[[p]], "R"), full.names = TRUE)) {
    ex <- parse(f, keep.source = TRUE)
    pd <- getParseData(ex)
    hits <- pd[pd$token == "SYMBOL_FUNCTION_CALL" &
                 pd$text %in% c("stopifnot", "match.arg"), ]
    sr <- attr(ex, "srcref")
    for (k in seq_len(nrow(hits))) {
      ln <- hits$line1[k]
      i <- which(vapply(sr, function(s) s[1] <= ln && s[3] >= ln, TRUE))
      top <- ex[[i]]
      owner <- if (is.call(top) && is.name(top[[2]])) as.character(top[[2]])
               else "?"
      call_txt <- getParseText(pd, pd$parent[pd$id == hits$parent[k]])
      rows[[length(rows) + 1L]] <- data.frame(
        pkg = p, file = basename(f), line = ln, fn = hits$text[k],
        owner = owner,
        exported = owner %in% exports || owner %in% ls(s3),
        call = gsub("\\s+", " ", call_txt), stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)
options(width = 200)
print(out, right = FALSE)
cat("\ntotal", nrow(out), " stopifnot", sum(out$fn == "stopifnot"),
    " match.arg", sum(out$fn == "match.arg"),
    " in an exported function or registered method", sum(out$exported), "\n")
