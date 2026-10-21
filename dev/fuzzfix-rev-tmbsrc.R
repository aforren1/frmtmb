# Reviewer: how TMB/RTMB build the integrate= transform (GK spec and
# the node count of marginal_gk).
for (ns in c("TMB", "RTMB")) {
  e <- asNamespace(ns)
  for (nm in ls(e, all.names = TRUE)) {
    f <- get(nm, e)
    if (!is.function(f)) next
    src <- deparse(f)
    if (any(grepl("marginal_gk|GK\\(|TransformADFunObject", src))) {
      cat("====", ns, "::", nm, "\n")
      cat(src[grepl("gk|GK|integrate|Transform|adaptive|which", src)],
          sep = "\n")
    }
  }
}
