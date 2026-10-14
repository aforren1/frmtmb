# Lane postfix: totals of one dev/postfix-run.sh tag from its per-file
# outputs: every file whose RESULT line shows a failure, an error, an
# escaped warning or a load error, the files with no RESULT line, the
# totals, and which core library the files read.
#   Rscript dev/postfix-summarize.R <tag> <job list>
a <- commandArgs(TRUE)
tag <- a[1]
jobs <- utils::read.table(a[2], col.names = c("p", "f"),
                          stringsAsFactors = FALSE)
out <- file.path("dev/postfix-tests", tag)
tot <- c(pass = 0, fail = 0, err = 0, skip = 0, warn = 0)
ran <- 0L
libs <- character()
for (i in seq_len(nrow(jobs))) {
  o <- file.path(out, paste0(jobs$p[i], "--",
                             sub("[.]R$", "", basename(jobs$f[i])), ".txt"))
  x <- if (file.exists(o)) readLines(o, warn = FALSE) else character()
  r <- utils::tail(grep("^RESULT", x, value = TRUE), 1L)
  libs <- c(libs, sub(".*; frmtmb [0-9.]+ ", "",
                      utils::head(grep("^lib:", x, value = TRUE), 1L)))
  if (!length(r)) {
    cat("NO RESULT", jobs$p[i], jobs$f[i], "\n")
    next
  }
  ran <- ran + 1L
  v <- vapply(names(tot), function(k) {
    as.numeric(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
  }, 0)
  if (grepl("LOADERROR", r)) cat(jobs$p[i], r, "\n")
  v[is.na(v)] <- 0
  tot <- tot + v
  if (any(v[c("fail", "err", "warn")] > 0)) cat(jobs$p[i], r, "\n")
}
cat(sprintf("TOTAL %s: files %d of %d; %s\n", tag, ran, nrow(jobs),
            paste(names(tot), tot, sep = "=", collapse = " ")))
cat("core library read:\n")
print(table(libs))
