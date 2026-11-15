# Route (b) of items 561/4846: does RTMB forward a Gauss-Kronrod spec's
# fields (TMB's gk_config: dx, ytol, ...) to TMBad's marginal_gk? Prints
# the RTMB code that reads `integrate`.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("RTMB")
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns)
  if (!is.function(f)) next
  src <- deparse(f)
  hit <- grep("integrate|marginal_gk|GK|gk_config|dx|ytol", src)
  if (length(hit)) {
    cat("=====", nm, "\n")
    lo <- max(1, min(hit) - 3); hi <- min(length(src), max(hit) + 3)
    cat(src[lo:hi], sep = "\n")
  }
}
cat("---- TMB::GK\n")
print(TMB::GK)
print(packageVersion("RTMB")); print(packageVersion("TMB"))
