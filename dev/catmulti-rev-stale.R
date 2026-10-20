# Reviewer: which functions of the worktree's R/ differ from the frmtmb
# installed in the lane library (installed 11:34; R/ edited 12:07).
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
cat("installed from:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-catmulti/R"
env <- new.env()
for (f in list.files(wt, full.names = TRUE, pattern = "[.]R$")) {
  exprs <- parse(f, keep.source = FALSE)
  for (e in exprs) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                       identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) {
      nm <- as.character(e[[2]])
      fn <- eval(e[[3]], baseenv())
      if (!exists(nm, envir = ns, inherits = FALSE)) {
        cat("MISSING in installed:", basename(f), nm, "\n")
        next
      }
      inst <- get(nm, envir = ns)
      if (!is.function(inst)) next
      a <- deparse(body(fn)); b <- deparse(body(inst))
      fa <- deparse(formals(fn)); fb <- deparse(formals(inst))
      if (!identical(a, b) || !identical(fa, fb)) {
        cat("DIFFERS:", basename(f), nm, "\n")
        d <- setdiff(a, b); d2 <- setdiff(b, a)
        cat("  source-only lines:\n"); cat(paste0("    ", head(d, 30)), sep = "\n")
        cat("  installed-only lines:\n"); cat(paste0("    ", head(d2, 30)), sep = "\n")
      }
    }
  }
}
cat("done\n")
