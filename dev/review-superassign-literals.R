# Review of lane superassign: function literals in core R/ that sit in a
# top-level object that is not itself a function (a list registry), and
# so are outside the guard's namespace scan.
cnt <- function(e) {
  if (!is.call(e)) return(0L)
  n <- if (identical(e[[1L]], as.name("function"))) 1L else 0L
  for (a in as.list(e)[-1L]) if (!missing(a)) n <- n + cnt(a)
  n
}
cat("lib: none (source scan of", normalizePath("R"), ")\n")
tot <- 0L
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  for (e in parse(f, keep.source = FALSE)) {
    if (!is.call(e) || !as.character(e[[1L]])[1L] %in% c("<-", "=")) next
    rhs <- e[[3L]]
    if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) next
    n <- cnt(rhs)
    if (n) {
      tot <- tot + n
      cat(sprintf("  %-22s %-28s %d\n", basename(f), deparse1(e[[2L]]), n))
    }
  }
}
cat("RESULT function literals outside a top-level function:", tot, "\n")
