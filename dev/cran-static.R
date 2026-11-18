# Lane cran: static checks CRAN's human review makes and R CMD check
# does not, on the worktree's sources. One R process; urlchecker reads
# the network (public URLs only, nothing posted).
#
#   Rscript dev/cran-static.R <worktree> > dev/cran-log/static.txt
local({
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- commandArgs(trailingOnly = TRUE)[1]
d <- read.dcf(file.path(root, "DESCRIPTION"))

cat("== Title case\n")
ttl <- unname(d[, "Title"])
cat("Title:      ", ttl, "\n")
cat("toTitleCase:", tools::toTitleCase(ttl), "\n")
cat("identical:  ", identical(ttl, tools::toTitleCase(ttl)), "\n")

cat("\n== NEWS.md parse (as CRAN's news() reads it)\n")
db <- tools:::.build_news_db_from_package_NEWS_md(
  file.path(root, "NEWS.md"))
cat("entries:", nrow(db), " versions:", length(unique(db$Version)),
    " bad:", sum(attr(db, "bad") %||% FALSE), "\n")
cat("NEWS.md bytes:", file.size(file.path(root, "NEWS.md")), "\n")

cat("\n== cat()/print() calls outside print/format/summary/show methods\n")
ns <- asNamespace("frmtmb")
fns <- Filter(function(n) is.function(get(n, ns)), ls(ns, all.names = TRUE))
hits <- character()
for (n in fns) {
  if (grepl("^(print|format|summary|show|toString)[.]", n)) next
  g <- codetools::findGlobals(get(n, ns), merge = FALSE)$functions
  u <- intersect(g, c("cat", "print", "writeLines"))
  if (length(u)) hits <- c(hits, paste0(n, ": ", paste(u, collapse = ",")))
}
cat(length(hits), "functions\n")
cat(hits, sep = "\n")

cat("\n== exported functions\n")
ex <- getNamespaceExports("frmtmb")
cat(length(ex), "exports\n")

cat("\n== urlchecker\n")
r <- tryCatch(urlchecker::url_check(root, progress = FALSE),
              error = function(e) {
                cat("url_check error:", conditionMessage(e), "\n")
                NULL
              })
if (!is.null(r)) {
  cat("urls with a problem:", nrow(r), "\n")
  if (nrow(r)) print(as.data.frame(r)[, c("URL", "Status", "Message",
                                          "New")], right = FALSE)
}
})
