# Reviewer of lane surface: library paths by arm.
#   "lane": the reviewer's install of the current worktree
#           (C:/Users/adf44/source/r/surface-rev-lib), then rellib-r6
#   "base": rellib-r6 alone (frmtmb 0.68.1)
rev_env <- function(arm = c("lane", "base")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/surface-rev-lib",
              "C:/Users/adf44/source/r/rellib-r6",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(
    R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
    FRMTMB_STAN_CACHE =
      "C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-rev-stan-cache",
    NOT_CRAN = "true")
  invisible(arm)
}

rev_show <- function(tag, expr) {
  w <- character(); m <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x))
    invokeRestart("muffleMessage")
  }),
  error = function(e) structure(conditionMessage(e), class = "rev_err"))
  cat(sprintf("%-46s %s\n", tag,
              if (inherits(r, "rev_err")) {
                paste("ERROR:", substr(gsub("\n", " ", r), 1, 150))
              } else paste("ok", paste(class(r), collapse = "/"))))
  for (x in unique(w)) cat("    W:", substr(gsub("\n", " ", x), 1, 150), "\n")
  for (x in unique(m)) cat("    M:", substr(gsub("\n", " ", x), 1, 150), "\n")
  invisible(r)
}
