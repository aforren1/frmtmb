.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib", "C:/Users/adf44/source/r/rellib-r11", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb.spatial")
fns <- Filter(is.function, mget(ls(ns, all.names = TRUE), ns))
g <- unique(unlist(lapply(fns, function(f) codetools::findGlobals(f, merge = FALSE)$functions)))
unres <- g[!vapply(g, function(x) exists(x, envir = ns), TRUE)]
cat("unresolved:", unres, "\n")
fx <- getNamespaceExports("frmtmb")
from_frmtmb <- g[vapply(g, function(x) exists(x, envir = ns) && identical(environmentName(environment(get(x, envir = ns))), "frmtmb"), TRUE)]
cat("functions from frmtmb namespace:", from_frmtmb, "\n")
cat("of which not exported:", setdiff(from_frmtmb, fx), "\n")
# inner closures: slots registered
reg <- frmtmb:::cs_entry("car")
inner <- unlist(lapply(Filter(is.function, c(reg, reg$term, reg$expand)), function(f) codetools::findGlobals(f, merge = FALSE)$functions))
cat("car entry closures call:", sort(unique(inner)), "\n")
cat("exports of frmtmb.spatial:", getNamespaceExports(ns), "\n")
