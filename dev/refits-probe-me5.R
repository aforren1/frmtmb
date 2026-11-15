.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("marginaleffects")
hits <- Filter(function(f) {
  fn <- get(f, ns)
  is.function(fn) &&
    any(grepl("not known to be supported", deparse(fn), fixed = TRUE))
}, ls(ns, all.names = TRUE))
print(hits)
for (h in hits) print(get(h, ns))
