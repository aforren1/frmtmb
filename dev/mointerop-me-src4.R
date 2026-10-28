.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in ls(ns, all.names = TRUE)) { g <- get(f, ns); if (is.function(g)) { b <- deparse(g); i <- grep("sanity_dots[(]|marginaleffects_safe", b); if (length(i)) { cat("=====", f, "\n"); cat(b[unique(pmax(1, rep(i, each = 5) + -2:2))], sep = "\n") } } }
