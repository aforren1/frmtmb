# Route (b) of 561/4846: where TMB reads an integrate spec's fields.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("TMB")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns)
  if (!is.function(f)) next
  src <- deparse(f)
  hit <- grep("integrate|marginal_gk|gk_config|ytol|\bdx\b", src)
  if (length(hit)) {
    cat("=====", nm, "\n")
    for (h in hit) cat(h, ":", src[max(1, h - 2):min(length(src), h + 6)], sep = "\n")
  }
}
