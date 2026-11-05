# Sum the RESULT lines of a suite directory: generated counts for the record.
d <- commandArgs(TRUE)[[1]]
fs <- Sys.glob(file.path(d, "*.txt"))
ln <- unlist(lapply(fs, function(f) {
  grep("^RESULT", readLines(f, warn = FALSE), value = TRUE)
}))
g <- function(k) {
  vapply(strsplit(ln, " ", fixed = TRUE), function(p) {
    as.integer(sub(paste0(k, "="), "", p[startsWith(p, paste0(k, "="))],
                   fixed = TRUE))
  }, 1L)
}
cat(sprintf("%s: files %d of %d; pass %d fail %d error %d skip %d warn %d\n",
            d, length(ln), length(fs), sum(g("pass")), sum(g("fail")),
            sum(g("error")), sum(g("skip")), sum(g("warn"))))
bad <- ln[g("fail") > 0 | g("error") > 0 | g("warn") > 0]
if (length(bad)) writeLines(substr(bad, 1, 100))
sk <- ln[g("skip") > 0]
if (length(sk)) {
  cat("with skips:\n")
  writeLines(substr(sk, 1, 100))
}
miss <- fs[!vapply(fs, function(f) {
  any(grepl("^RESULT", readLines(f, warn = FALSE)))
}, NA)]
if (length(miss)) {
  cat("no RESULT:\n")
  writeLines(miss)
}
