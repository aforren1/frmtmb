# Fix round 1: is the installed library the tree under review? Every top-level function
# of R/ against the installed namespace, by deparse.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
cmp_pkg <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  diffs <- character(); n <- 0L
  for (f in list.files(dir, pattern = "[.]R$", full.names = TRUE)) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (!is.call(e) || !(identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) || !is.name(e[[2]])) next
      if (!is.call(e[[3]]) || !identical(e[[3]][[1]], as.name("function"))) next
      nm <- as.character(e[[2]])
      if (!exists(nm, envir = ns, inherits = FALSE)) { diffs <- c(diffs, paste("missing", nm)); next }
      g <- get(nm, envir = ns)
      if (!is.function(g)) next
      src <- eval(e[[3]], envir = baseenv())
      n <- n + 1L
      if (!identical(deparse(body(src)), deparse(body(utils::removeSource(g))))) diffs <- c(diffs, paste(basename(f), nm))
    }
  }
  cat(sprintf("RB2 %s lib %s functions compared %d differing %d\n", pkg, dirname(system.file(package = pkg)), n, length(diffs)))
  print(diffs)
}
cmp_pkg("frmtmb", file.path(REV_WT, "R"))
cmp_pkg("frmtmb.sample", file.path(REV_WT, "extensions/frmtmb.sample/R"))
