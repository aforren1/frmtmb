# Lane code8 final verification: does the installed frmtmb in
# wt-code8-lib carry the worktree's R/fit.R and R/se-check.R?
.libPaths(c("C:/Users/adf44/source/r/wt-code8-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-code8"
ns <- asNamespace("frmtmb")
cat("installed:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
for (file in c("R/fit.R", "R/se-check.R")) {
  e <- new.env()
  exprs <- parse(file.path(wt, file), keep.source = FALSE)
  for (x in exprs) {
    if (is.call(x) && (identical(x[[1]], as.name("<-")) ||
                         identical(x[[1]], as.name("="))) &&
          is.name(x[[2]]) && is.call(x[[3]]) &&
          identical(x[[3]][[1]], as.name("function"))) {
      assign(as.character(x[[2]]), eval(x[[3]], e), envir = e)
    }
  }
  nm <- ls(e, all.names = TRUE)
  diff <- character()
  for (n in nm) {
    if (!exists(n, envir = ns, inherits = FALSE)) {
      diff <- c(diff, paste0(n, " (missing)"))
      next
    }
    a <- deparse(get(n, envir = e))
    b <- deparse(removeSource(get(n, envir = ns)))
    if (!identical(a, b)) diff <- c(diff, n)
  }
  cat(file, ": ", length(nm), " functions; differ: ",
      if (length(diff)) paste(diff, collapse = ", ") else "none", "\n",
      sep = "")
}
for (n in c("headroom_ritz_flat", "newton_room", "headroom_off_flat")) {
  cat(n, "in namespace:", exists(n, envir = ns, inherits = FALSE), "\n")
}
