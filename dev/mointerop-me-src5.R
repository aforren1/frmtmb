.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
for (f in c("get_predict.default", "get_predict.glmmTMB", "get_predict.merMod", "pad_newdata", "complete_levels")) { cat("=====", f, "\n"); if (exists(f, ns)) print(get(f, ns)) }
for (f in ls(ns, all.names = TRUE)) { g <- get(f, ns); if (is.function(g)) { b <- deparse(g); i <- grep("rowid < 0|rowid > 0|rowid_dedup|padding", b); if (length(i)) { cat("=====", f, "\n"); cat(b[unique(pmax(1, rep(i, each = 5) + -2:2))], sep = "\n") } } }
