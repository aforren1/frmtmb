# Reviewer of lane loglik: is the lane library built from the worktree
# as it stands? Every function of the changed R files, source against
# installed. Output: dev/review-loglik-log/sync.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-loglik"
chk <- function(pkg, files) {
  ns <- asNamespace(pkg)
  for (f in files) {
    env <- new.env()
    sys.source(file.path(wt, f), envir = env, keep.source = FALSE)
    nd <- 0L
    for (nm in ls(env, all.names = TRUE)) {
      a <- get(nm, env)
      if (!is.function(a) || !exists(nm, ns, inherits = FALSE)) next
      b <- get(nm, ns)
      if (!identical(deparse(body(a)), deparse(body(b))) ||
          !identical(deparse(formals(a)), deparse(formals(b)))) {
        cat("DIFFERS", pkg, f, nm, "\n")
        nd <- nd + 1L
      }
    }
    cat(pkg, f, "differences:", nd, "\n")
  }
}
chk("frmtmb", c("R/autocor-loglik.R", "R/autocor-seam.R"))
chk("frmtmb.sample", "extensions/frmtmb.sample/R/loo.R")
chk("frmtmb.spatial", "extensions/frmtmb.spatial/R/sar.R")
