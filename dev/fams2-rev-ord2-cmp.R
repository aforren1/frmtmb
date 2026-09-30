od <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out"
b <- readRDS(file.path(od, "ord2-base.rds")); l <- readRDS(file.path(od, "ord2-lane.rds"))
ns <- 0; nd <- 0
walk <- function(x, y, path) {
  if (!is.null(x$value) || !is.null(x$warnings) || identical(names(x), c("value", "warnings"))) {
    same <- identical(x, y)
    isErr <- is.character(x$value) && length(x$value) == 1 && startsWith(x$value, "ERROR")
    cat(sprintf("%-4s %-40s %s\n", if (same) "same" else "DIFF", path,
                if (isErr) substr(x$value, 1, 90) else ""))
    if (same) ns <<- ns + 1 else { nd <<- nd + 1; print(head(all.equal(x, y), 5)) }
    return(invisible())
  }
  for (k in names(x)) walk(x[[k]], y[[k]], paste(path, k, sep = "/"))
}
walk(b, l, "")
cat(sprintf("identical %d, differing %d\n", ns, nd))
