## Reviewer: does the worker's installed build match the worktree source?
## Compares the deparsed bodies of every function the diff touches
## against a fresh parse of the worktree R/ sources.
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("libPaths:\n"); print(.libPaths())
cat("frmtmb from:", find.package("frmtmb"), "\n")
cat("version:", as.character(packageVersion("frmtmb")), "\n")
suppressMessages(library(frmtmb))

WT <- "C:/Users/adf44/source/r/frmtmb-wt-gradcheck"
env <- new.env()
for (f in list.files(file.path(WT, "R"), pattern = "[.]R$",
                     full.names = TRUE)) {
  try(sys.source(f, envir = env, keep.source = FALSE), silent = TRUE)
}
fns <- c("fit_outer_box", "grad_bound_active", "grad_headroom",
         "grad_verdict", "grad_warning_msg", "check_convergence",
         "diagnose", "fit_assembled", "autoscale_prefit")
for (nm in fns) {
  a <- tryCatch(get(nm, envir = asNamespace("frmtmb")),
                error = function(e) NULL)
  b <- tryCatch(get(nm, envir = env), error = function(e) NULL)
  if (is.null(a)) { cat(sprintf("%-20s ABSENT in build\n", nm)); next }
  if (is.null(b)) { cat(sprintf("%-20s ABSENT in worktree\n", nm)); next }
  same <- identical(deparse(a), deparse(b))
  cat(sprintf("%-20s body identical: %s\n", nm, same))
}
cat("--- he() availability probe ---\n")
cat("sessionInfo RTMB:", as.character(packageVersion("RTMB")), "\n")
