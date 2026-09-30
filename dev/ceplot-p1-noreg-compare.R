# Lane ceplot punch 1: a copy of the reviewer's comparison: compare the two arms of dev/ceplot-rev-noreg.R
# with identical(), naming what differs where they do.
#   Rscript dev/ceplot-rev-noreg-compare.R > dev/ceplot-rev-log/noreg-compare.txt
lg <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-log/p1-"
a <- readRDS(paste0(lg, "noreg-base.rds"))
b <- readRDS(paste0(lg, "noreg-lane.rds"))
nsame <- 0L
ndiff <- 0L
nerr <- 0L
for (m in union(names(a), names(b))) {
  for (k in union(names(a[[m]]), names(b[[m]]))) {
    x <- a[[m]][[k]]
    y <- b[[m]][[k]]
    if (is.character(x) && length(x) == 1L && startsWith(x, "ERROR")) {
      nerr <- nerr + 1L
    }
    if (identical(x, y)) {
      nsame <- nsame + 1L
      if (is.character(x) && length(x) == 1L) {
        cat(sprintf("%-10s %-13s identical (both: %s)\n", m, k,
                    substr(x, 1, 90)))
      }
      next
    }
    ndiff <- ndiff + 1L
    why <- if (is.character(x) || is.character(y)) {
      paste0("base: ", substr(paste(x), 1, 100), " | lane: ",
             substr(paste(y), 1, 100))
    } else if (is.list(x) && is.list(y)) {
      sub <- vapply(union(names(x), names(y)), function(n) {
        if (identical(unclass(x[[n]]), unclass(y[[n]]))) {
          if (identical(attributes(x[[n]]), attributes(y[[n]]))) "" else
            paste0(n, "[attributes: ",
                   paste(setdiff(union(names(attributes(x[[n]])),
                                       names(attributes(y[[n]]))),
                                 names(Filter(isTRUE, Map(identical,
                                   attributes(x[[n]])[intersect(names(attributes(x[[n]])), names(attributes(y[[n]])))],
                                   attributes(y[[n]])[intersect(names(attributes(x[[n]])), names(attributes(y[[n]])))])))),
                         collapse = ","), "]")
        } else paste0(n, "[values: max diff ",
                      tryCatch(format(max(abs(as.numeric(unlist(x[[n]])) -
                                                as.numeric(unlist(y[[n]]))),
                                          na.rm = TRUE), digits = 3),
                               error = function(e) "?"), "]")
      }, "")
      top <- if (!identical(attributes(x), attributes(y))) " + top attributes" else ""
      paste(c(sub[nzchar(sub)], top), collapse = " ")
    } else {
      paste("max diff", format(max(abs(x - y), na.rm = TRUE), digits = 3))
    }
    cat(sprintf("%-10s %-13s DIFFERS: %s\n", m, k, why))
  }
}
cat(sprintf("identical %d, differ %d (errors on base %d)\n", nsame, ndiff,
            nerr))
