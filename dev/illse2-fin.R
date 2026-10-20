# Lane illse2, resumed: which loaded namespace registers a finalizer
# `function(ptr)`, the one whose "promise already under evaluation"
# error dev/illse2-p4wide.R 150 printed after the library restore.
#   Rscript dev/illse2-fin.R
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (p in loadedNamespaces()) {
  ns <- asNamespace(p)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get0(nm, envir = ns, inherits = FALSE)
    if (!is.function(f)) next
    src <- tryCatch(paste(deparse(body(f)), collapse = "\n"),
                    error = function(e) "")
    if (grepl("reg.finalizer", src, fixed = TRUE)) {
      cat(p, "::", nm, "\n")
      cat(grep("reg.finalizer", strsplit(src, "\n")[[1]], value = TRUE),
          sep = "\n")
    }
  }
}
for (p in c("RTMB", "TMB", "Matrix", "frmtmb")) {
  cat(p, format(packageVersion(p)), find.package(p),
      packageDescription(p)$Packaged %||% "",
      packageDescription(p)$Built, "\n")
}
