# Reviewer: does the lane build match the worktree source, function by
# function? A stale build would make every other measurement here a
# measurement of the wrong code.
LIB <- "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"

src_env <- function(dir) {
  e <- new.env(parent = globalenv())
  fs <- sort(list.files(dir, pattern = "[.][Rr]$", full.names = TRUE))
  for (f in fs) {
    tryCatch(sys.source(f, envir = e, keep.source = FALSE),
             error = function(c) message("source err ", basename(f), ": ",
                                         conditionMessage(c)))
  }
  e
}

cmp_pkg <- function(pkg, dir) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ns <- asNamespace(pkg)
  e <- src_env(dir)
  nms <- ls(e, all.names = TRUE)
  diffs <- character(0); missing <- character(0)
  for (nm in nms) {
    o <- get(nm, envir = e)
    if (!is.function(o)) next
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      missing <- c(missing, nm); next
    }
    p <- get(nm, envir = ns)
    if (!is.function(p)) { diffs <- c(diffs, paste0(nm, " [not fn]")); next }
    a <- paste(deparse(o), collapse = "\n")
    b <- paste(deparse(p), collapse = "\n")
    if (!identical(a, b)) diffs <- c(diffs, nm)
  }
  # names in the namespace not in source (generated / installed extras)
  extra <- setdiff(
    Filter(function(n) is.function(get(n, envir = ns)),
           ls(ns, all.names = TRUE)),
    nms)
  cat("== ", pkg, " ==\n", sep = "")
  cat("  source functions:", sum(vapply(nms, function(n)
    is.function(get(n, envir = e)), TRUE)), "\n")
  cat("  differing bodies:", length(diffs), "\n")
  if (length(diffs)) cat("    ", paste(diffs, collapse = ", "), "\n")
  cat("  in source, not in namespace:", length(missing), "\n")
  if (length(missing)) cat("    ", paste(missing, collapse = ", "), "\n")
  cat("  in namespace, not in source:", length(extra), "\n")
  if (length(extra) && length(extra) < 40)
    cat("    ", paste(extra, collapse = ", "), "\n")
  invisible(NULL)
}

cmp_pkg("frmtmb", file.path(TREE, "R"))
cmp_pkg("frmtmb.sample",
        file.path(TREE, "extensions/frmtmb.sample/R"))
cat("\ninstalled versions\n")
for (p in c("frmtmb", "frmtmb.sample", "frmtmb.spline", "frmtmb.eam",
            "frmtmb.coupling", "frmtmb.latent", "frmtmb.learn",
            "frmtmb.ode"))
  cat(sprintf("  %-18s %s  built %s\n", p,
              as.character(packageVersion(p)),
              utils::packageDescription(p)$Built))
