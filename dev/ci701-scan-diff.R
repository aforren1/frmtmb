# Lane ci701: compare scan runs (dev/ci701-scan.sh) file by file: the
# RESULT line, every DETAIL line (non-passing expectation) and every
# TEST pass count. Prints what differs between the runs.
#
#   Rscript dev/ci701-scan-diff.R <run> <run> [<run> ...]
runs <- commandArgs(TRUE)
logd <- "C:/Users/adf44/source/r/frmtmb-wt-ci701/dev/ci701-log"
read_run <- function(r) {
  fs <- list.files(file.path(logd, paste0(r, "-suite")), pattern = "[.]txt$",
                   full.names = TRUE)
  out <- list()
  for (f in fs) {
    l <- readLines(f, warn = FALSE)
    keep <- grep("^(RESULT|DETAIL|TEST)", l, value = TRUE)
    keep <- gsub("[0-9]+[.][0-9]+(e-?[0-9]+)?", "#", keep)
    out[[basename(f)]] <- keep
  }
  out
}
d <- lapply(runs, read_run)
files <- sort(unique(unlist(lapply(d, names))))
ndiff <- 0L
for (f in files) {
  v <- lapply(d, function(x) x[[f]])
  if (any(vapply(v, is.null, NA))) {
    cat("MISSING in some run:", f, "\n"); ndiff <- ndiff + 1L; next
  }
  if (!all(vapply(v[-1], identical, NA, v[[1]]))) {
    ndiff <- ndiff + 1L
    cat("DIFF", f, "\n")
    for (i in seq_along(runs)) {
      only <- setdiff(v[[i]], Reduce(intersect, v))
      for (x in only) cat("  ", runs[i], ":", substr(x, 1, 200), "\n")
    }
  }
}
cat("files compared", length(files), "; differing", ndiff, "\n")
