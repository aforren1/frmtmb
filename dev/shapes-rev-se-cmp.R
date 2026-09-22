# Compare the two captures. A difference is reported as the maximum
# relative change, never as a printed rounding, because %.6f renders
# 1.8e-11 as zero.
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
b <- readRDS(file.path(TREE, "dev/shapes-rev-se-base.rds"))
l <- readRDS(file.path(TREE, "dev/shapes-rev-se-lane.rds"))

reldiff <- function(x, y) {
  if (inherits(x, "revErr") && inherits(y, "revErr"))
    return(list(kind = "both-error", v = NA_real_,
                note = paste0("base: ", x$msg, " | lane: ", y$msg)))
  if (inherits(x, "revErr"))
    return(list(kind = "base-error-only", v = NA_real_, note = x$msg))
  if (inherits(y, "revErr"))
    return(list(kind = "LANE-ERROR-ONLY", v = NA_real_, note = y$msg))
  if (is.list(x) && !is.data.frame(x)) {
    return(list(kind = "list",
                v = max(vapply(names(x), function(k)
                  reldiff(x[[k]], y[[k]])$v %||% NA_real_, 0), na.rm = TRUE)))
  }
  if (is.data.frame(x)) { x <- as.matrix(x[sapply(x, is.numeric)])
                          y <- as.matrix(y[sapply(y, is.numeric)]) }
  if (is.character(x))
    return(list(kind = if (identical(x, y)) "same-chr" else "CHR-DIFF",
                v = 0, note = paste(setdiff(union(x, y),
                                            intersect(x, y)),
                                    collapse = ",")))
  if (!identical(dim(x), dim(y)))
    return(list(kind = "DIM-DIFF", v = NA_real_,
                note = paste0(paste(dim(x), collapse = "x"), " -> ",
                              paste(dim(y), collapse = "x"))))
  nmx <- dimnames(x); nmy <- dimnames(y)
  namediff <- !identical(nmx, nmy)
  d <- abs(as.numeric(x) - as.numeric(y))
  s <- pmax(abs(as.numeric(x)), 1e-12)
  fin <- is.finite(d)
  v <- if (any(fin)) max((d / s)[fin]) else NA_real_
  nab <- sum(!is.finite(as.numeric(x))) ; nal <- sum(!is.finite(as.numeric(y)))
  list(kind = if (namediff) "names-differ" else "same-names", v = v,
       note = sprintf("nonfinite base=%d lane=%d", nab, nal))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat(sprintf("%-10s %-14s %-18s %12s  %s\n",
            "fit", "quantity", "kind", "max rel diff", "note"))
for (fnm in names(b)) {
  if (fnm == "importance") {
    r <- reldiff(b[[fnm]], l[[fnm]])
    cat(sprintf("%-10s %-14s %-18s %12.3e  %s\n", "-", fnm, r$kind,
                r$v, r$note %||% ""))
    next
  }
  for (q in names(b[[fnm]])) {
    r <- reldiff(b[[fnm]][[q]], l[[fnm]][[q]])
    flag <- if (!is.na(r$v) && r$v > 1e-8) "  <== MOVED" else ""
    if (grepl("DIFF|ERROR", r$kind)) flag <- "  <== STRUCT"
    cat(sprintf("%-10s %-14s %-18s %12.3e  %s%s\n", fnm, q, r$kind,
                r$v, r$note %||% "", flag))
  }
}
