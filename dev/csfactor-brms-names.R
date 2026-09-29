# What does brms 2.23.0 call a cs() factor's coefficients in
# variables()? Read the renaming rule out of brms's own namespace rather
# than sampling a model for it.
#   Rscript dev/csfactor-brms-names.R > dev/csfactor-log/brms-names.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
ns <- asNamespace("brms")
fns <- ls(ns, all.names = TRUE)
hit <- fns[vapply(fns, function(f) {
  ob <- get(f, envir = ns)
  is.function(ob) && any(grepl("bcs", deparse(ob), fixed = TRUE))
}, TRUE)]
cat("functions mentioning bcs:", paste(hit, collapse = ", "), "\n")
for (f in hit) {
  src <- deparse(get(f, envir = ns))
  ln <- grep("bcs", src, fixed = TRUE)
  cat("\n== ", f, " ==\n", sep = "")
  cat(paste(src[sort(unique(c(ln, pmax(1, ln - 3), pmin(length(src),
                                                        ln + 3))))],
            collapse = "\n"), "\n")
}
