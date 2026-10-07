# Reviewer: dump every function of a namespace as deparsed text, to
# compare the build the lane's suite ran on with the current sources.
# usage: Rscript surface-rev-fndump.R <lib> <pkg> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace(a[2])
nm <- ls(ns, all.names = TRUE)
out <- lapply(nm, function(n) {
  f <- get(n, ns)
  if (is.function(f)) paste(deparse(f, control = NULL), collapse = "\n")
  else paste(deparse(f), collapse = "\n")
})
names(out) <- nm
saveRDS(out, a[3])
cat(a[2], "from", find.package(a[2]), length(out), "objects\n")
