# Lane illse2, resumed: closures with a single formal `ptr`, in every
# loaded namespace and in the environments RTMB keeps.
#   Rscript dev/illse2-fin2.R
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
walk <- function(x, depth = 0) {
  if (depth > 6) return(character())
  out <- character()
  if (is.function(x) && identical(names(formals(x)), "ptr")) {
    out <- paste(deparse(x), collapse = " ")
  }
  if (is.function(x) && !is.primitive(x)) {
    b <- body(x)
    out <- c(out, unlist(lapply(as.list(b), function(e) {
      if (is.call(e)) walk_call(e, depth + 1) else character()
    })))
  }
  out
}
walk_call <- function(e, depth) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("function[(]ptr[)]", s)) substr(s, 1, 300) else character()
}
for (p in loadedNamespaces()) {
  ns <- asNamespace(p)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get0(nm, envir = ns, inherits = FALSE)
    if (!is.function(f)) next
    r <- unique(walk(f))
    if (length(r)) cat(p, "::", nm, ":", r[1], "\n\n")
  }
}
