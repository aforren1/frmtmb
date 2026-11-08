# Review of lane mfxse, second short check: is the installed build the source?
# R/interop.R (14:10:23) is newer than the library (14:10:18). Compare every
# function R/interop.R and R/autocor.R define with the installed namespace.
.libPaths(c("C:/Users/adf44/source/r/wt-mfxse-lib", "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
for (f in c("interop.R", "autocor.R")) {
  ex <- parse(file.path("C:/Users/adf44/source/r/frmtmb-wt-mfxse/R", f), keep.source = FALSE)
  env <- new.env(parent = ns)
  nd <- 0; n <- 0
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]])) {
      nm <- as.character(e[[2]])
      v <- tryCatch(eval(e[[3]], env), error = function(err) NULL)
      if (!is.function(v) || !exists(nm, ns, inherits = FALSE)) next
      n <- n + 1
      a <- deparse(v); b <- deparse(get(nm, ns))
      if (!identical(a, b)) { nd <- nd + 1; cat("  differs:", nm, "\n") }
    }
  }
  cat(f, ": functions compared", n, ", differing", nd, "\n")
}
