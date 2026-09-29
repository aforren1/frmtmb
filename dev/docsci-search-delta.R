# The exact difference between the committed root search index and the
# one the final build produced, and the same for the sitemap.
one <- function(x) if (is.null(x) || !length(x)) "" else as.character(x)[1]
a <- jsonlite::fromJSON("docs/search.json", simplifyVector = FALSE)
b <- jsonlite::fromJSON("dev/docsci-site/search.json",
                        simplifyVector = FALSE)
pa <- vapply(a, function(x) one(x$path), character(1))
pb <- vapply(b, function(x) one(x$path), character(1))
cat("search entries: committed ", length(a), ", new ", length(b), "\n",
    sep = "")
cat("subsite paths: committed ",
    sum(grepl("/frmtmb/frmtmb[.]", pa)), ", new ",
    sum(grepl("/frmtmb/frmtmb[.]", pb)), "\n", sep = "")
cat("only in committed:\n")
cat(paste0("  ", setdiff(pa, pb)), sep = "\n")
cat("\nonly in new:\n")
cat(paste0("  ", setdiff(pb, pa)), sep = "\n")

loc <- function(p) {
  x <- readLines(p, warn = FALSE)
  m <- unlist(regmatches(x, gregexpr("<loc>[^<]+</loc>", x)))
  sub("</loc>$", "", sub("^<loc>", "", m))
}
la <- loc("docs/sitemap.xml")
lb <- loc("dev/docsci-site/sitemap.xml")
cat("\nsitemap locs: committed ", length(la), ", new ", length(lb), "\n",
    sep = "")
cat("subsite locs: committed ", sum(grepl("/frmtmb/frmtmb[.]", la)),
    ", new ", sum(grepl("/frmtmb/frmtmb[.]", lb)), "\n", sep = "")
only <- setdiff(la, lb)
cat("only in committed: ", length(only), "\n", sep = "")
cat("of those, under the CORE reference: ",
    sum(grepl("/frmtmb/reference/", only)), "\n", sep = "")
cat("first 5:\n")
cat(paste0("  ", head(only, 5)), sep = "\n")
cat("\nonly in new: ", length(setdiff(lb, la)), "\n", sep = "")
if (length(setdiff(lb, la))) {
  cat(paste0("  ", head(setdiff(lb, la), 10)), sep = "\n")
}
