# Reviewer, claim 1: where the stored frames differ, list the leaf
# paths; and list which elements errored (on both arms alike).
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/"
a <- readRDS(paste0(wt, "aterms2-rev-01-base.rds"))
b <- readRDS(paste0(wt, "aterms2-rev-01-lane.rds"))
leafdiff <- function(x, y, path = "") {
  if (identical(x, y)) return(character(0))
  if (is.list(x) && is.list(y)) {
    nx <- names(x)
    ny <- names(y)
    if (is.null(nx) || is.null(ny)) {
      if (length(x) != length(y)) return(paste0(path, " [length ",
                                                length(x), " vs ",
                                                length(y), "]"))
      return(unlist(lapply(seq_along(x), function(i)
        leafdiff(x[[i]], y[[i]], paste0(path, "[[", i, "]]")))))
    }
    out <- character(0)
    for (k in setdiff(union(nx, ny), "")) {
      if (!k %in% nx) out <- c(out, paste0(path, "$", k, " only lane: ",
                                           paste(format(y[[k]]),
                                                 collapse = " ")))
      else if (!k %in% ny) out <- c(out, paste0(path, "$", k, " only base"))
      else out <- c(out, leafdiff(x[[k]], y[[k]], paste0(path, "$", k)))
    }
    return(out)
  }
  paste0(path, " base=", paste(format(x), collapse = " "), " lane=",
         paste(format(y), collapse = " "))
}
for (m in names(a)) {
  d <- leafdiff(a[[m]]$frame, b[[m]]$frame, "frame")
  if (length(d)) cat(m, ":\n  ", paste(substr(d, 1, 200),
                                        collapse = "\n   "), "\n")
}
cat("\nerrored elements per model (identical on both arms):\n")
for (m in names(a)) {
  e <- names(a[[m]])[vapply(a[[m]], function(x)
    is.character(x) && length(x) == 1L && startsWith(x, "ERROR"), TRUE)]
  cat(sprintf("%-24s %s\n", m, paste(e, collapse = ",")))
}
cat("\nsample error messages:\n")
msgs <- unlist(lapply(a, function(r) Filter(function(x)
  is.character(x) && length(x) == 1L && startsWith(x, "ERROR"), r)))
print(head(sort(table(substr(msgs, 1, 110)), decreasing = TRUE), 25))
