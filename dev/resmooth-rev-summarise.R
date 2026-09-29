# Reviewer: emit the suite counts from the RESULT files, refusing to
# print totals while any launched file has no RESULT line.
dirs <- c("dev/resmooth-rev-suite", "dev/resmooth-rev-base")
for (dd in dirs) {
  if (!dir.exists(dd)) next
  want <- readLines(file.path(dd, "files.lst"), warn = FALSE)
  want <- basename(want[nzchar(want)])
  got <- list.files(dd, pattern = "\\.R\\.txt$")
  missing <- setdiff(paste0(want, ".txt"), got)
  rows <- list()
  noresult <- character(0)
  for (g in got) {
    ln <- grep("^RESULT", readLines(file.path(dd, g), warn = FALSE),
               value = TRUE)
    if (!length(ln)) { noresult <- c(noresult, g); next }
    n <- function(k) as.integer(sub(paste0(".*\\b", k, "=(\\d+).*"), "\\1",
                                    ln[1L]))
    rows[[g]] <- c(pass = n("pass"), fail = n("fail"), error = n("error"),
                   skip = n("skip"), warn = n("warn"))
  }
  cat("==", dd, "==\n")
  cat("  launched:", length(want), "| with a RESULT line:", length(rows),
      "\n")
  if (length(missing)) cat("  NEVER RAN:", paste(missing, collapse = " "),
                           "\n")
  if (length(noresult)) cat("  NO RESULT LINE:",
                            paste(noresult, collapse = " "), "\n")
  if (length(missing) || length(noresult)) {
    cat("  totals withheld until every file reports\n\n"); next
  }
  M <- do.call(rbind, rows)
  cat(sprintf("  pass=%d fail=%d error=%d skip=%d warn=%d\n",
              sum(M[, "pass"]), sum(M[, "fail"]), sum(M[, "error"]),
              sum(M[, "skip"]), sum(M[, "warn"])))
  bad <- rownames(M)[M[, "fail"] > 0 | M[, "error"] > 0]
  cat("  files with a fail or error:", length(bad),
      if (length(bad)) paste0("(", paste(bad, collapse = ", "), ")") else "",
      "\n")
  sk <- rownames(M)[M[, "skip"] > 0]
  if (length(sk)) {
    for (s in sk) cat("    skip:", s, M[s, "skip"], "\n")
  }
  cat("\n")
}
