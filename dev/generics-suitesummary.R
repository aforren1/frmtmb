# Summarise a one-file-per-process suite run into a block that cannot
# be mistaken for a clean sweep.
#
# The counts come from the LOGS, not from the launcher, because a
# launcher tally counts files it started rather than files that
# finished.  A log with no FILE line is reported as NO RESULT rather
# than skipped silently.
#
#   Rscript dev/generics-suitesummary.R <outdir> [label]
av <- commandArgs(trailingOnly = TRUE)
dir <- av[1]
lab <- if (length(av) >= 2L) av[2] else dir
logs <- sort(list.files(dir, pattern = "[.]log$", full.names = TRUE))

num <- function(ln, k) {
  v <- grep(paste0("^", k, " "), ln, value = TRUE)
  if (!length(v)) return(NA_integer_)
  as.integer(trimws(sub(paste0("^", k), "", v[1])))
}

rows <- lapply(logs, function(p) {
  ln <- readLines(p, warn = FALSE)
  data.frame(file = basename(p),
             blocks = num(ln, "BLOCKS"), assert = num(ln, "ASSERT"),
             pass = num(ln, "PASS"), fail = num(ln, "FAIL"),
             error = num(ln, "ERROR"), skip = num(ln, "SKIP"),
             warn = num(ln, "WARN"), secs = NA_real_,
             stringsAsFactors = FALSE)
})
d <- do.call(rbind, rows)
noresult <- d$file[is.na(d$assert)]
ok <- d[!is.na(d$assert), , drop = FALSE]

cat("```\n")
cat(sprintf("== suite: %s ==\n", lab))
cat(sprintf("test files with a result   %d of %d\n", nrow(ok), nrow(d)))
cat(sprintf("test blocks                %d\n", sum(ok$blocks)))
cat(sprintf("assertions                 %d\n", sum(ok$assert)))
cat(sprintf("pass                       %d\n", sum(ok$pass)))
cat(sprintf("FAIL                       %d\n", sum(ok$fail)))
cat(sprintf("ERROR                      %d\n", sum(ok$error)))
cat(sprintf("skip                       %d\n", sum(ok$skip)))
cat(sprintf("warning                    %d\n", sum(ok$warn)))
if (length(noresult)) {
  cat("\nNO RESULT (the process wrote no count):\n")
  for (f in noresult) cat("  ", f, "\n")
}
bad <- ok[ok$fail > 0 | ok$error > 0, , drop = FALSE]
if (nrow(bad)) {
  cat("\nfiles with a failure or an error:\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %-40s fail %d error %d\n", bad$file[i], bad$fail[i],
                bad$error[i]))
  }
} else {
  cat("\nno file reported a failure or an error\n")
}
cat("```\n")
