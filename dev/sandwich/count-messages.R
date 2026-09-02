# srr G5.2a: every condition message raised in R/ is unique. Parse R/,
# collect the literal string fragments of every stop() call, and require
# no duplicates. Reproduces 529 at v0.35.0 and 554 with the
# cluster-robust batch.
#
# Run from the package root: Rscript dev/sandwich/count-messages.R
msgs <- character(0)
walk <- function(e) {
  if (is.call(e)) {
    f <- e[[1L]]
    if (is.name(f) && as.character(f) == "stop") {
      lits <- unlist(lapply(as.list(e)[-1L], function(a) {
        if (is.character(a)) a else NULL
      }))
      if (length(lits)) msgs <<- c(msgs, paste(lits, collapse = ""))
    }
  }
  if (is.call(e) || is.pairlist(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) {
      if (!is.null(e[[i]])) try(walk(e[[i]]), silent = TRUE)
    }
  }
  invisible(NULL)
}
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  for (e in parse(f)) walk(e)
}
cat("stop() calls carrying literal text:", length(msgs), "\n")
cat("distinct messages:", length(unique(msgs)), "\n")
dups <- unique(msgs[duplicated(msgs)])
if (length(dups)) {
  cat("DUPLICATES:\n")
  writeLines(paste0("  ", dups))
  quit(status = 1L)
}
