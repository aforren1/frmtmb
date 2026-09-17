# Reviewer: does the lane library hold the worktree's R sources?
# Every top-level function assignment in R/*.R is evaluated alone and
# compared with the installed namespace binding, body and formals.
lib <- "C:/Users/adf44/source/r/priorform-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
cat("frmtmb from", find.package("frmtmb"), "\n")
files <- list.files("C:/Users/adf44/source/r/frmtmb-wt-priorform/R",
                    "[.]R$", full.names = TRUE)
n_ok <- 0L; bad <- character(0); missing <- character(0)
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (!is.call(e) || !(identical(e[[1]], as.name("<-")) ||
                         identical(e[[1]], as.name("=")))) next
    if (!is.name(e[[2]]) && !is.character(e[[2]])) next
    rhs <- e[[3]]
    if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
    nm <- as.character(e[[2]])
    fn <- eval(rhs, envir = baseenv())
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      missing <- c(missing, nm); next
    }
    got <- get(nm, envir = ns, inherits = FALSE)
    if (!is.function(got)) { missing <- c(missing, nm); next }
    same <- identical(deparse(body(fn)), deparse(body(got))) &&
      identical(deparse(formals(fn)), deparse(formals(got)))
    if (same) n_ok <- n_ok + 1L else bad <- c(bad, paste0(basename(f), ":", nm))
  }
}
cat("functions identical:", n_ok, "\n")
cat("differ:", length(bad), paste(bad, collapse = " "), "\n")
cat("not a function in namespace (redefined or active):", length(missing),
    paste(head(missing, 20), collapse = " "), "\n")
