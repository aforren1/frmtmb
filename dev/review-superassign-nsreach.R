# Review of lane superassign: how much package code the guard's
# namespace scan cannot reach. It reads only the namespace's top-level
# functions; a function literal held in a top-level list or environment
# (a registry such as autocor_core, covstruct_core) or in an S4 method
# table is not walked. Under R CMD check on a tarball the source scan
# skips, so there the namespace scan is the only guard.
.libPaths(c("C:/Users/adf44/source/r/wt-superassign-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkgs <- c("frmtmb", "frmtmb.ode", "frmtmb.sample", "frmtmb.survival")
for (p in pkgs) suppressMessages(loadNamespace(p))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
for (p in pkgs) {
  ns <- asNamespace(p)
  own <- topenv(ns)
  seen <- new.env()
  found <- character(0)
  visit <- function(x, path, depth) {
    if (depth > 8L) return()
    if (is.function(x) && !is.primitive(x)) {
      if (identical(topenv(environment(x)), own)) found <<- c(found, path)
      return()
    }
    if (is.environment(x)) {
      if (isNamespace(x) || identical(x, globalenv()) ||
            identical(x, emptyenv()) || identical(x, baseenv())) return()
      key <- format(x)
      if (!is.null(seen[[key]])) return()
      seen[[key]] <- TRUE
      for (nm in ls(x, all.names = TRUE)) {
        v <- tryCatch(get(nm, envir = x, inherits = FALSE),
                      error = function(e) NULL)
        visit(v, paste0(path, "$", nm), depth + 1L)
      }
      return()
    }
    if (is.list(x)) {
      nms <- names(x) %||% rep("", length(x))
      for (i in seq_along(x)) {
        visit(x[[i]], paste0(path, "[[", if (nzchar(nms[i])) nms[i] else i,
                             "]]"), depth + 1L)
      }
    }
  }
  top_fun <- 0L
  for (nm in ls(ns, all.names = TRUE)) {
    v <- tryCatch(get(nm, envir = ns, inherits = FALSE),
                  error = function(e) NULL)
    if (is.function(v)) {
      top_fun <- top_fun + 1L
      next
    }
    visit(v, nm, 0L)
  }
  s4 <- ls(ns, all.names = TRUE, pattern = "^\\.__T__")
  holders <- unique(sub("[$[].*", "", found))
  cat(sprintf("%-16s top-level functions %d; functions held in lists or environments %d (in %d objects); S4 method tables %d\n",
              p, top_fun, length(found), length(holders), length(s4)))
  if (length(holders)) cat("   holders:", paste(head(holders, 40), collapse = ", "), "\n")
}
