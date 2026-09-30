# Lane ceplot: the counts dev/ceplot-findings.md quotes, generated from
# the logs rather than typed.
#   Rscript dev/ceplot-findings-gen.R <tag> [<tag> ...]
# Each tag is a dev/ceplot-log/<tag>/ directory of per-file logs
# written by dev/ceplot-par.sh. Prints, per tag, the number of files,
# the RESULT lines that are not clean, the totals, and the files with
# no RESULT line.
tags <- commandArgs(trailingOnly = TRUE)
root <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-log"
num <- function(l, k) {
  as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", l))
}
for (tag in tags) {
  fs <- list.files(file.path(root, tag), full.names = TRUE)
  res <- vapply(fs, function(f) {
    l <- grep("^RESULT", readLines(f, warn = FALSE), value = TRUE)
    if (length(l)) sub("\r$", "", l[length(l)]) else NA_character_
  }, "")
  cat("## ", tag, ": ", length(fs), " files, ", sum(!is.na(res)),
      " with a RESULT line\n", sep = "")
  ok <- res[!is.na(res)]
  load_err <- grepl("LOADERROR", ok)
  ok2 <- ok[!load_err]
  for (arm in unique(sub("^RESULT ([a-z]+) .*", "\\1", ok2))) {
    for (p in unique(sub("^RESULT [a-z]+ ([^ ]+) .*", "\\1", ok2))) {
      sel <- ok2[startsWith(ok2, paste0("RESULT ", arm, " ", p, " "))]
      if (!length(sel)) next
      cat(sprintf(
        "%s %s: %d files, expectations %d, failed %d, error %d, skipped %d, warning %d, passed %d\n",
        arm, p, length(sel), sum(num(sel, "tests")), sum(num(sel, "failed")),
        sum(num(sel, "error")), sum(num(sel, "skipped")),
        sum(num(sel, "warning")), sum(num(sel, "passed"))))
    }
  }
  bad <- ok2[num(ok2, "failed") > 0 | num(ok2, "error") > 0 |
               num(ok2, "warning") > 0]
  if (length(bad)) {
    cat("not clean:\n")
    for (i in seq_along(bad)) {
      cat("  ", bad[i], "\n", sep = "")
      f <- names(bad)[i]
      for (l in grep("^  (fails|warns):", readLines(f, warn = FALSE),
                     value = TRUE)) {
        cat("    ", sub("\r$", "", trimws(l)), "\n", sep = "")
      }
    }
  }
  if (any(load_err)) cat("LOADERROR:\n", paste0("  ", ok[load_err], "\n"))
  if (any(is.na(res))) {
    cat("no RESULT line:\n", paste0("  ", basename(fs[is.na(res)]), "\n"))
  }
  cat("\n")
}
