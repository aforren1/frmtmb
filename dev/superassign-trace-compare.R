# Compare the base and lane traces written by dev/superassign-trace.R.
# Prints, per test file, the RESULT line of each arm, the number of
# traced calls, and every call whose function or hash differs; then the
# calls per traced function summed over all files (lane arm).
#
# Usage: Rscript dev/superassign-trace-compare.R <trace-out dir>

dir <- commandArgs(trailingOnly = TRUE)[1L]
lanes <- sort(list.files(dir, pattern = "--lane[.]tsv$", full.names = TRUE))
total <- list()
n_diff_files <- 0L
for (lf in lanes) {
  bf <- sub("--lane[.]tsv$", "--base.tsv", lf)
  key <- sub("--lane[.]tsv$", "", basename(lf))
  res <- function(arm) {
    txt <- readLines(file.path(dir, paste0(key, "--", arm, ".txt")),
                     warn = FALSE)
    r <- grep("^RESULT", txt, value = TRUE)
    if (length(r)) sub("^RESULT [^ ]+ ", "", r[length(r)]) else "NO RESULT"
  }
  if (!file.exists(bf)) {
    cat(sprintf("%s: base trace missing\n", key))
    next
  }
  L <- utils::read.delim(lf, colClasses = "character")
  B <- utils::read.delim(bf, colClasses = "character")
  same_n <- nrow(L) == nrow(B)
  m <- min(nrow(L), nrow(B))
  dif <- which(L$fn[seq_len(m)] != B$fn[seq_len(m)] |
                 L$hash[seq_len(m)] != B$hash[seq_len(m)])
  if (!same_n || length(dif)) n_diff_files <- n_diff_files + 1L
  cat(sprintf("%s: calls base=%d lane=%d differing=%d%s | base %s | lane %s\n",
              key, nrow(B), nrow(L), length(dif),
              if (same_n) "" else " LENGTH", res("base"), res("lane")))
  for (i in utils::head(dif, 10L)) {
    cat(sprintf("    call %d: base %s %s, lane %s %s\n", i, B$fn[i],
                B$hash[i], L$fn[i], L$hash[i]))
  }
  tb <- table(L$fn)
  for (nm in names(tb)) total[[nm]] <- (total[[nm]] %||% 0L) + tb[[nm]]
}
cat(sprintf("FILES %d, with any difference %d\n", length(lanes),
            n_diff_files))
cat("CALLS PER FUNCTION (lane arm, all files)\n")
for (nm in sort(names(total))) cat(sprintf("  %-26s %d\n", nm, total[[nm]]))
