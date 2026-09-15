# Lane eamhier: read the per-replicate records into one data frame.
#
# One file per replicate, one tab-separated key=value line in each. The
# reader is shared by every summary script so that no two of them can
# disagree about what a record means.

read_records <- function(dirs) {
  fs <- unlist(lapply(dirs, function(p) {
    list.files(p, pattern = "[.]tsv$", full.names = TRUE)
  }))
  rows <- lapply(fs, function(f) {
    ln <- readLines(f, warn = FALSE)
    ln <- ln[nzchar(ln)]
    if (!length(ln)) return(NULL)
    kv <- strsplit(ln[[1L]], "\t", fixed = TRUE)[[1L]]
    nm <- sub("=.*$", "", kv)
    vl <- sub("^[^=]*=", "", kv)
    out <- as.list(trimws(vl))
    names(out) <- nm
    out$file <- f
    out
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  nms <- unique(unlist(lapply(rows, names)))
  d <- do.call(rbind, lapply(rows, function(r) {
    r <- r[nms]
    names(r) <- nms
    as.data.frame(lapply(r, function(x) if (is.null(x)) NA else x),
                  stringsAsFactors = FALSE)
  }))
  # everything numeric except the fields that are words. A record whose
  # `status` is not "ok" carries a message, which is why `msg` is here.
  for (nm in setdiff(nms, c("arm", "bound", "status", "unbounded",
                            "vc_names", "vc_sd", "msg", "file",
                            "pdHess", "fit_sv"))) {
    d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  }
  d
}

# Wilson, because the normal interval on a proportion near 1 runs past
# it and 60 replicates is not large.
wilson <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  c(((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}
