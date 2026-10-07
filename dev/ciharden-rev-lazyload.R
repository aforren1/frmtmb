# Reviewer: does base_r_function() load a namespace when the name IS a
# function of an unloaded package? Forcing a promise from a lazy-load
# database unserializes the closure, whose environment is the package
# namespace, and R finds a namespace reference by loading it.
#   Rscript dev/ciharden-rev-lazyload.R <lib|base> <name>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.sample))
probe <- c("tcltk", "stats4", "splines", "tools", "compiler", "grid",
           "parallel")
before <- vapply(probe, isNamespaceLoaded, NA)
ans <- frmtmb.sample:::base_r_function(a[2])
after <- vapply(probe, isNamespaceLoaded, NA)
cat(sprintf("%-6s name=%-12s answer=%s newly loaded: %s\n",
            if (a[1] == "base") "base" else "lane", a[2], ans,
            paste(probe[after & !before], collapse = ",")))
# the same question asked the obvious way, for the answer it should give
truth <- any(vapply(probe, function(p) {
  exists(a[2], envir = asNamespace(p), inherits = FALSE) &&
    is.function(get(a[2], envir = asNamespace(p)))
}, NA))
cat("        truth over the probe packages:", truth, "\n")
