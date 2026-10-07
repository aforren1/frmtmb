# Reviewer: does the installed lane build match the worktree source?
# Compares every top-level function of the named files, by deparse.
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden"
chk <- function(pkg, files) {
  ns <- asNamespace(pkg)
  cat(pkg, "from", dirname(find.package(pkg)), "\n")
  for (f in files) {
    ex <- parse(file.path(wt, f), keep.source = FALSE)
    nd <- 0L
    nf <- 0L
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
            is.name(e[[2]])) {
        nm <- as.character(e[[2]])
        v <- tryCatch(eval(e[[3]], envir = new.env(parent = ns)),
                      error = function(err) NULL)
        if (!is.function(v) || !exists(nm, ns, inherits = FALSE)) next
        g <- get(nm, ns)
        if (!is.function(g)) next
        nf <- nf + 1L
        if (!identical(deparse(body(v)), deparse(body(g))) ||
              !identical(deparse(formals(v)), deparse(formals(g)))) {
          nd <- nd + 1L
          cat("  DIFFERS:", f, nm, "\n")
        }
      }
    }
    cat(" ", f, ":", nf, "functions,", nd, "differ\n")
  }
}
chk("frmtmb", c("R/frame.R", "R/predict.R", "R/utils.R"))
chk("frmtmb.sample", c("extensions/frmtmb.sample/R/sample.R",
                       "extensions/frmtmb.sample/R/methods-draws.R"))
