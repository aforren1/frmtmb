## Does the lane library hold the worktree's R/ code? Compares every
## function defined in R/*.R (parsed from source) with the installed
## namespace's copy of it.
.libPaths(c("C:/Users/adf44/source/r/wt-mocode-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
lib <- commandArgs(trailingOnly = TRUE)[1L]
if (!is.na(lib)) .libPaths(c(lib, .libPaths()))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(dirname(system.file(package = "frmtmb"))), "\n")
src <- "C:/Users/adf44/source/r/frmtmb-wt-mocode/R"
diffs <- character(0)
for (f in list.files(src, full.names = TRUE, pattern = "[.]R$")) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        is.name(e[[2L]]) && is.call(e[[3L]]) &&
        identical(e[[3L]][[1L]], as.name("function"))) {
      nm <- as.character(e[[2L]])
      if (!exists(nm, envir = ns, inherits = FALSE)) next
      g <- get(nm, envir = ns)
      if (!is.function(g)) next
      a <- deparse(body(eval(e[[3L]])))
      b <- deparse(body(g))
      if (!identical(a, b)) diffs <- c(diffs, paste(basename(f), nm))
    }
  }
}
cat("functions differing:", length(diffs), "\n")
cat(diffs, sep = "\n")
