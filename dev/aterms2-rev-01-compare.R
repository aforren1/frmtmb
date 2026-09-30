# Reviewer, claim 1: compare the base and lane outputs of
# dev/aterms2-rev-01-regress.R with identical(), element by element.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/"
a <- readRDS(paste0(wt, "aterms2-rev-01-base.rds"))
b <- readRDS(paste0(wt, "aterms2-rev-01-lane.rds"))
stopifnot(identical(names(a), names(b)))
nel <- 0L
ndiff <- 0L
nerr <- 0L
for (m in names(a)) {
  els <- union(names(a[[m]]), names(b[[m]]))
  diffs <- character(0)
  for (e in els) {
    nel <- nel + 1L
    x <- a[[m]][[e]]
    y <- b[[m]][[e]]
    if (is.character(x) && length(x) == 1L && startsWith(x, "ERROR")) {
      nerr <- nerr + 1L
    }
    if (!identical(x, y)) {
      ndiff <- ndiff + 1L
      diffs <- c(diffs, e)
    }
  }
  errs <- vapply(els, function(e) {
    x <- a[[m]][[e]]
    is.character(x) && length(x) == 1L && startsWith(x, "ERROR")
  }, TRUE)
  cat(sprintf("%-24s elements %3d  differ %2d  base-errors %2d %s\n", m,
              length(els), length(diffs), sum(errs),
              if (length(diffs)) paste(diffs, collapse = ",") else ""))
  for (e in diffs) {
    cat("   ", e, ": base ", substr(paste(format(a[[m]][[e]]),
                                          collapse = " "), 1, 160), "\n")
    cat("   ", e, ": lane ", substr(paste(format(b[[m]][[e]]),
                                          collapse = " "), 1, 160), "\n")
  }
}
cat("TOTAL models", length(a), "elements", nel, "differ", ndiff,
    "error-elements (same on both unless listed)", nerr, "\n")
