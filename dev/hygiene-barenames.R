# Item 769: list the names a test file reads that a session with the
# package attached cannot see, but the package namespace (or a
# namespace it imports from) can. Those are the names that make the
# file runnable only under load_all() or test_check().
#
# Usage: Rscript dev/hygiene-barenames.R <arm: mine|base> <pkg> <file>...
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]; p <- a[2]; files <- a[-(1:2)]
.libPaths(c(if (arm == "mine") "C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(testthat)
  library(p, character.only = TRUE)
})
ns <- asNamespace(p)
imp <- parent.env(ns)
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  nm <- unique(all.names(ex))
  # names the file or its helpers assign are local, not internals
  helpers <- list.files(dirname(f), "^helper.*[.][Rr]$", full.names = TRUE)
  local_def <- unique(unlist(lapply(c(f, helpers), function(h) {
    e <- parse(h, keep.source = FALSE)
    out <- character(0)
    walk <- function(x) {
      if (missing(x)) return(invisible())
      if (is.call(x)) {
        fn <- x[[1L]]
        if (is.symbol(fn) && as.character(fn) %in%
            c("<-", "=", "<<-", "assign") && length(x) >= 2L) {
          lhs <- x[[2L]]
          if (is.symbol(lhs)) out <<- c(out, as.character(lhs))
          if (is.character(lhs)) out <<- c(out, lhs)
        }
        if (is.symbol(fn) && as.character(fn) == "function") {
          out <<- c(out, names(x[[2L]]))
        }
        args <- as.list(x)[-1L]
        keep <- !vapply(args, function(a) {
          is.symbol(a) && !nzchar(as.character(a))
        }, logical(1))
        lapply(args[keep], walk)
      } else if (is.pairlist(x)) {
        out <<- c(out, names(x))
      }
    }
    for (i in seq_along(e)) walk(e[[i]])
    out
  })))
  hidden <- nm[!vapply(nm, exists, logical(1), envir = globalenv()) &
               (vapply(nm, exists, logical(1), envir = ns, inherits = FALSE) |
                vapply(nm, exists, logical(1), envir = imp,
                       inherits = FALSE))]
  hidden <- setdiff(hidden, local_def)
  where <- vapply(hidden, function(s) {
    if (exists(s, envir = ns, inherits = FALSE)) p else
      environmentName(environment(get(s, envir = imp)))
  }, character(1))
  cat("FILE ", f, " hidden=", length(hidden), "\n", sep = "")
  for (s in hidden) cat("  ", s, " (", where[[s]], ")\n", sep = "")
}
