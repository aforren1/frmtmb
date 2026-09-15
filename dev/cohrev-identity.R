## Load-independent check: is the LANE's installed frmtmb.coupling
## doing the same work as the BASE build's? If every object in the
## namespace deparses identically, the ?cross_wishart example cannot
## have got slower because of this lane.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1L]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.coupling))
cat("path:", dirname(dirname(getNamespaceInfo("frmtmb.coupling",
                                              "path"))), "\n")
cat("where:", getNamespaceInfo("frmtmb.coupling", "path"), "\n")
ns <- asNamespace("frmtmb.coupling")
nms <- sort(ls(ns, all.names = TRUE))
h <- vapply(nms, function(n) {
  o <- get(n, envir = ns)
  paste(deparse(o, control = c("useSource", "keepInteger",
                               "showAttributes")), collapse = "\n")
}, character(1))
out <- args[2L]
saveRDS(list(nms = nms, h = h), out)
cat("objects:", length(nms), "\n")
