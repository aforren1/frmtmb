.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in ls(ns, all.names = TRUE)) { g <- get(f, ns); if (is.function(g)) { b <- paste(deparse(g), collapse = "\n"); if (grepl("known to be supported", b)) { cat("=====", f, "\n"); print(g) } } }
