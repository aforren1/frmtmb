# Short review of lane alapse: is the lane library built from the
# worktree's source? Every top-level function of the changed R files,
# deparsed, against the installed namespace.
.libPaths(c("C:/Users/adf44/source/r/wt-alapse-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
files <- c("R/ald-fixed.R", "R/fit.R", "R/confint.R", "R/se-check.R",
           "R/ald-laplace.R", "R/families.R", "R/compat.R", "R/interop.R")
for (f in files) {
  e <- new.env()
  sys.source(f, envir = e, keep.source = FALSE)
  nms <- ls(e, all.names = TRUE)
  diff <- character()
  for (n in nms) {
    a <- get(n, e)
    if (!exists(n, ns, inherits = FALSE)) { diff <- c(diff, paste0(n, "(missing)")); next }
    b <- get(n, ns)
    if (is.function(a)) {
      same <- identical(deparse(body(a)), deparse(body(b))) &&
        identical(deparse(formals(a)), deparse(formals(b)))
    } else same <- isTRUE(all.equal(a, b))
    if (!same) diff <- c(diff, n)
  }
  cat(sprintf("%-18s objects %4d differ %d %s\n", f, length(nms),
              length(diff), paste(head(diff, 10), collapse = " ")))
}
