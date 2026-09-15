# Which namespace pulls coda (and any other owner) in when only
# frmtmb.sample is loaded; an owner shadow must not break a dependency.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setHook(packageEvent("coda", "onLoad"), function(...) {
  cat("coda loaded; call stack namespaces:\n")
  print(sys.calls()[grepl("loadNamespace", vapply(sys.calls(),
    function(x) paste(deparse(x)[1], collapse = ""), ""))])
})
suppressMessages(loadNamespace("frmtmb.sample"))
own <- c("posterior", "loo", "rstantools", "bayesplot", "coda",
         "bridgesampling", "gratia", "brms")
cat("owners loaded after loadNamespace(frmtmb.sample):",
    intersect(own, loadedNamespaces()), "\n")
