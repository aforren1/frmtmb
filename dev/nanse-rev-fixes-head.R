# Reviewer: lane fixes' flat-case table (dev/fixes-p2-table.R in the
# fixes worktree, body copied below this header by
# dev/nanse-rev-fixes-build.sh), run on the trial merge: per model, every
# warning frm() raises, classified, so a fit warned twice shows up.
#   Rscript dev/nanse-rev-fixes-table.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
kind <- function(m) {
  if (grepl("^Standard errors are not available", m)) return("SE")
  if (grepl("are not identified: at the optimum", m)) return("NLFLAT")
  if (grepl("^Optimizer did not report|gradient", m)) return("CONV")
  if (grepl("^Some standard errors are not finite", m)) return("VCOV")
  if (grepl("^Hessian is not positive", m)) return("PDHESS")
  "OTHER"
}
run <- function(lab, expr, truth) {
  w <- character()
  r <- tryCatch({
    fit <- suppressMessages(withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }))
    wv <- character()
    se <- withCallingHandlers(fixef(fit)[, "Est.Error"],
      warning = function(x) {
        wv <<- c(wv, conditionMessage(x)); invokeRestart("muffleWarning")
      })
    k <- vapply(w, kind, "")
    kv <- vapply(wv, kind, "")
    sprintf("%-4s conv %d  frm(): %-14s fixef(): %-6s finite SE %d/%d",
            truth, fit$opt$convergence,
            if (length(k)) paste(k, collapse = "+") else "none",
            if (length(kv)) paste(kv, collapse = "+") else "none",
            sum(is.finite(se)), length(se))
  }, error = function(e) paste(truth, "ERROR:",
                               substr(conditionMessage(e), 1, 100)))
  cat(sprintf("%-46s %s\n", lab, r))
}
