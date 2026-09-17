# Lane wt-conditions: the shapes of stop()/warning()/message() calls,
# so every non-plain site is decided by hand before the rewrite.
roots <- c("R", file.path(list.dirs("extensions", recursive = FALSE), "R"))
kinds <- c("stop", "warning", "message", "warningCondition",
           "errorCondition", "stopifnot")
argnames <- list()
odd <- list()
nsq <- 0L
walk <- function(e, f) {
  if (!is.call(e)) return(invisible())
  h <- e[[1L]]
  if (is.call(h) && identical(h[[1L]], as.name("::")) &&
      as.character(h[[3L]]) %in% kinds) {
    nsq <<- nsq + 1L
    cat("NAMESPACED:", f, deparse(e)[1], "\n")
  }
  if (is.name(h) && as.character(h) %in% kinds) {
    k <- as.character(h)
    nms <- names(e); if (is.null(nms)) nms <- rep("", length(e))
    for (n in nms[-1]) if (nzchar(n))
      argnames[[paste(k, n)]] <<- c(argnames[[paste(k, n)]], f)
    un <- which(!nzchar(nms[-1])) + 1L
    first <- if (length(un)) e[[un[1]]] else NULL
    plain <- is.character(first) ||
      (is.call(first) && is.name(first[[1L]]) &&
         as.character(first[[1L]]) %in% c("paste0", "paste", "sprintf",
                                          "gettextf", "sQuote", "dQuote",
                                          "format", "toupper"))
    if (k %in% c("stop", "warning", "message") &&
        (!plain || length(un) == 0L)) {
      odd[[length(odd) + 1L]] <<- paste0(f, ": ", paste(deparse(e),
                                         collapse = " "))
    }
  }
  for (i in seq_along(e)) {
    a <- tryCatch(if (is.call(e[[i]])) e[[i]], error = function(err) NULL)
    if (!is.null(a)) walk(a, f)
  }
}
for (r in roots) for (f in list.files(r, "[.][Rr]$", full.names = TRUE))
  for (e in parse(f, keep.source = FALSE)) walk(e, f)
cat("namespaced:", nsq, "\n")
print(vapply(argnames, length, 1L))
for (n in names(argnames)) {
  if (!grepl("call[.]$", n)) cat(n, ":", unique(argnames[[n]]), "\n")
}
cat(length(odd), "non-literal first argument sites\n")
writeLines(substr(unlist(odd), 1, 220))
