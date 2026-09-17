## Recheck, priority 2: which installed interop packages read a family's
## link with `[[`, getElement or .subset2 (which the read-time `$` does
## not serve)? Scans deparsed function bodies of each namespace.
.libPaths(c("C:/Users/adf44/source/r/pinlib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("insight", "emmeans", "marginaleffects", "broom", "ggeffects", "modelbased", "performance", "parameters")) {
  if (!requireNamespace(p, quietly = TRUE)) { cat(p, ": not installed\n"); next }
  ns <- asNamespace(p); hits <- character()
  for (nm in ls(ns, all.names = TRUE)) {
    fn <- get(nm, ns); if (!is.function(fn)) next
    txt <- paste(deparse(fn), collapse = "\n")
    if (any(vapply(c("[[\"link\"]]", "[[\"linkinv\"]]", "[[\"linkfun\"]]", "getElement(family", ".subset2(family"), function(s) grepl(s, txt, fixed = TRUE), TRUE))) hits <- c(hits, nm)
  }
  cat(p, ": functions reading [[\"link\"]]-style:", if (length(hits)) paste(hits, collapse = ", ") else "none", "\n")
}
# positive control: the same scan must find `$link` readers
ns <- asNamespace("insight"); n <- 0L
for (nm in ls(ns, all.names = TRUE)) { fn <- get(nm, ns); if (is.function(fn) && grepl("$link", paste(deparse(fn), collapse = ""), fixed = TRUE)) n <- n + 1L }
cat("control: insight functions containing `$link`:", n, "\n")
