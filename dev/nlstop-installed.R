# Lane nlstop: the installed core's functions equal the tree's R/ files.
.libPaths(c("C:/Users/adf44/source/r/wt-nlstop-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
e <- new.env(parent = ns)
for (f in c("R/fit-restart.R", "R/fit.R", "R/par-template.R", "R/families.R")) sys.source(f, e, keep.source = FALSE)
nm <- ls(e)
d <- nm[!vapply(nm, function(n) {
  a <- get(n, e); b <- get0(n, ns, inherits = FALSE)
  if (is.function(a)) is.function(b) && identical(deparse(a), deparse(b)) else identical(a, b)
}, NA)]
cat("objects", length(nm), "differ", length(d), if (length(d)) paste(d, collapse = " "), "\n")
