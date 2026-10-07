# Reviewer of lane setier, final check: does the installed lane library
# hold the tree's R/se-check.R and R/fit.R? Every function defined in
# them is compared, deparsed, with the namespace's.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (f in c("R/se-check.R", "R/fit.R")) {
  e <- new.env()
  sys.source(f, envir = e, keep.source = FALSE)
  nm <- ls(e, all.names = TRUE)
  bad <- nm[!vapply(nm, function(n) {
    exists(n, ns, inherits = FALSE) &&
      identical(deparse(get(n, e)), deparse(get(n, ns)))
  }, NA)]
  cat(f, ": objects", length(nm), "differing", length(bad),
      paste(head(bad), collapse = " "), "\n")
}
