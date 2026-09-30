# Reviewer: where brms 2.23.0 reads its `offset` argument.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
for (f in ls(ns)) {
  x <- get(f, ns)
  if (!is.function(x)) next
  s <- deparse(x)
  i <- grep("offset", s)
  i <- i[grepl("offset = |!offset|isFALSE[(]offset|[(]offset[)]|use_offset",
               s[i])]
  if (length(i)) cat("==", f, "\n", paste(s[i], collapse = "\n"), "\n")
}
