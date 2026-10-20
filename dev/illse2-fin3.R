# Lane illse2, resumed: finalizers in TMB and RTMB's R code.
#   Rscript dev/illse2-fin3.R
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("TMB", "RTMB")) {
  ns <- asNamespace(p)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get0(nm, envir = ns, inherits = FALSE)
    if (!is.function(f)) next
    s <- strsplit(paste(deparse(f), collapse = "\n"), "\n")[[1]]
    i <- grep("finaliz|function[(]ptr", s)
    if (length(i)) {
      cat("==", p, "::", nm, "\n")
      cat(s[unique(pmin(length(s), c(i, i + 1)))], sep = "\n")
    }
  }
}
