.libPaths(c("C:/Users/adf44/source/r/pinlib",
             "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
fs <- sub("^[.]family_", "", ls(asNamespace("brms"), all.names = TRUE, pattern = "^[.]family_"))
for (f in fs) {
  g <- get(paste0(".family_", f), asNamespace("brms"))
  if (length(formals(g))) { cat(f, "TAKES ARGS\n"); next }
  i <- g()
  cat(sprintf("%-32s %-5s %s\n", f, i$type %||% "", paste(i$links, collapse=",")))
}
print(brms:::print.brmsfamily)
print(brms:::use_real); print(brms:::use_int); print(brms:::no_mixture)
