# Which frmtmb packages each extension names in its DESCRIPTION, by
# field, so the CI install order can be checked against them.
for (p in list.files("extensions", full.names = TRUE)) {
  d <- read.dcf(file.path(p, "DESCRIPTION"))
  for (f in intersect(c("Depends", "Imports", "Suggests"), colnames(d))) {
    x <- trimws(sub("[(].*", "", unlist(strsplit(d[, f], ","))))
    x <- x[startsWith(x, "frmtmb")]
    if (length(x)) cat(basename(p), f, ":", x, "\n")
  }
  if ("Additional_repositories" %in% colnames(d)) {
    cat(basename(p), "Additional_repositories :",
        d[, "Additional_repositories"], "\n")
  }
}
