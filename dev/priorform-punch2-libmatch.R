# Does the lane library hold the worktree's R sources, for core and for
# frmtmb.sample? A stale install would make every pin and suite below
# measure the wrong code.
#   Rscript dev/priorform-punch2-libmatch.R
lib <- "C:/Users/adf44/source/r/priorform-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
check_pkg <- function(pkg, rdir) {
  ns <- asNamespace(pkg)
  files <- list.files(rdir, "[.]R$", full.names = TRUE)
  n_ok <- 0L; bad <- character(0); missing <- character(0)
  for (f in files) {
    for (e in parse(f, keep.source = FALSE)) {
      if (!is.call(e) || !(identical(e[[1]], as.name("<-")) ||
                           identical(e[[1]], as.name("=")))) next
      if (!is.name(e[[2]]) && !is.character(e[[2]])) next
      rhs <- e[[3]]
      if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
      nm <- as.character(e[[2]])
      fn <- eval(rhs, envir = baseenv())
      got <- if (exists(nm, envir = ns, inherits = FALSE)) {
        get(nm, envir = ns, inherits = FALSE)
      }
      if (!is.function(got)) { missing <- c(missing, nm); next }
      same <- identical(deparse(body(fn)), deparse(body(got))) &&
        identical(deparse(formals(fn)), deparse(formals(got)))
      if (same) n_ok <- n_ok + 1L else {
        bad <- c(bad, paste0(basename(f), ":", nm))
      }
    }
  }
  cat(sprintf("LIBMATCH %s from %s: identical %d differ %d not-a-function %d\n",
              pkg, find.package(pkg), n_ok, length(bad), length(missing)))
  if (length(bad)) cat("  differ:", bad, "\n")
  if (length(missing)) cat("  not a function:", head(missing, 30), "\n")
  invisible(n_ok)
}
n1 <- check_pkg("frmtmb", file.path(root, "R"))
n2 <- check_pkg("frmtmb.sample",
                file.path(root, "extensions", "frmtmb.sample", "R"))
# the comparison must have compared something, or it proves nothing
stopifnot(n1 > 500L, n2 > 50L)
