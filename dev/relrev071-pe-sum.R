# Summary of dev/relrev071-pe.ps1 runs against the release record's
# reference run (dev/ciharden-log/r071-ref-suite/): per configuration,
# files with a RESULT line, any fail, error, skip or warning, and files
# whose pass count differs from the reference run's.
#   Rscript dev/relrev071-pe-sum.R <pe dir>
d <- commandArgs(trailingOnly = TRUE)[1]
refd <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/ciharden-log/r071-ref-suite"
rd <- function(f) {
  L <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", L, value = TRUE)
  if (!length(r)) return(NULL)
  r <- r[length(r)]
  v <- as.integer(sub(".*=", "", strsplit(r, " ")[[1]][3:7]))
  names(v) <- c("pass", "fail", "err", "skip", "warn")
  v
}
for (cfg in list.dirs(d, recursive = FALSE, full.names = FALSE)) {
  fs <- list.files(file.path(d, cfg), "\\.txt$")
  tot <- c(pass = 0, fail = 0, err = 0, skip = 0, warn = 0)
  bad <- character(); moved <- character(); nores <- character()
  for (f in fs) {
    v <- rd(file.path(d, cfg, f))
    if (is.null(v)) { nores <- c(nores, f); next }
    tot <- tot + v
    if (any(v[-1] > 0)) bad <- c(bad, sprintf("%s %s", f,
      paste(names(v), v, sep = "=", collapse = " ")))
    r <- if (file.exists(file.path(refd, f))) rd(file.path(refd, f))
    if (!is.null(r) && r[["pass"]] != v[["pass"]]) {
      moved <- c(moved, sprintf("%s pass %d (ref run %d)", f, v[["pass"]],
                                r[["pass"]]))
    }
  }
  cat(sprintf("== %s: %d files, %d with RESULT; pass %d fail %d err %d skip %d warn %d\n",
              cfg, length(fs), length(fs) - length(nores), tot[["pass"]],
              tot[["fail"]], tot[["err"]], tot[["skip"]], tot[["warn"]]))
  if (length(nores)) cat("  no RESULT:", paste(nores, collapse = ", "), "\n")
  for (b in bad) cat("  NONZERO", b, "\n")
  for (m in moved) cat("  PASS COUNT", m, "\n")
}
