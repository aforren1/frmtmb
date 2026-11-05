# Compare the two arms of dev/review-nlfix-a8.R: exact, then in ulps.
rd <- function(f) {
  l <- readLines(f)
  sp <- strsplit(l, "\t", fixed = TRUE)
  setNames(vapply(sp, function(x) paste(x[-(1:2)], collapse = "\t"), ""),
           vapply(sp, function(x) paste(x[1:2], collapse = "/"), ""))
}
a <- rd("dev/review-nlfix-a8-before.tsv")
b <- rd("dev/nlfix-a8-after.tsv")
cat("keys only before:", setdiff(names(a), names(b)), "\n")
cat("keys only after:", setdiff(names(b), names(a)), "\n")
k <- intersect(names(a), names(b))
same <- a[k] == b[k]
cat("identical:", sum(same), "of", length(k), "\n")
ncases <- length(unique(sub("/.*", "", k)))
cat("cases:", ncases, "\n")
for (nm in k[!same]) {
  va <- suppressWarnings(as.numeric(strsplit(a[[nm]], ",", fixed = TRUE)[[1]]))
  vb <- suppressWarnings(as.numeric(strsplit(b[[nm]], ",", fixed = TRUE)[[1]]))
  if (length(va) == length(vb) && !anyNA(va) && !anyNA(vb)) {
    cat(sprintf("%-28s max abs %.3g  max rel %.3g\n", nm, max(abs(va - vb)), max(abs(va - vb) / pmax(1e-300, abs(va)))))
  } else {
    cat(nm, "\n   before:", substr(a[[nm]], 1, 400), "\n   after: ", substr(b[[nm]], 1, 400), "\n")
  }
}
