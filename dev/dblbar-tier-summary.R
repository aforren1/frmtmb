# Lane dblbar: per-package totals of a tier log written by
# dev/dblbar-run-par.sh, and every file that is not clean, generated so
# that no count is typed. Also checks that every per-file log names the
# lane's core library.
# Usage: Rscript dev/dblbar-tier-summary.R <tier> [<core lib dir>]
a <- commandArgs(TRUE)
tier <- a[1]
want_lib <- if (length(a) > 1L) a[2] else ""
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
lg <- readLines(file.path(wt, "dev/dblbar-tests", paste0(tier, ".log")))
pkg <- NA_character_
rows <- list()
for (l in lg) {
  m <- regmatches(l, regexec("^== (\\S+) [0-9]+ files ==$", l))[[1]]
  if (length(m)) pkg <- m[2]
  if (startsWith(l, "RESULT ")) {
    f <- strsplit(l, " ")[[1]][2]
    s <- sub("^RESULT \\S+ ", "", l)
    n <- as.integer(regmatches(s, gregexpr("[0-9]+", s))[[1]])
    rows[[length(rows) + 1L]] <- data.frame(pkg = pkg, file = f,
                                            pass = n[1], fail = n[2],
                                            err = n[3], skip = n[4],
                                            warn = n[5])
  }
}
d <- do.call(rbind, rows)
tot <- aggregate(cbind(files = 1, pass, fail, err, skip, warn) ~ pkg, d, sum)
print(tot, row.names = FALSE)
cat("all:", nrow(d), "files; pass", sum(d$pass), "fail", sum(d$fail),
    "err", sum(d$err), "skip", sum(d$skip), "warn", sum(d$warn), "\n")
bad <- d[d$fail > 0 | d$err > 0 | d$warn > 0, ]
if (nrow(bad)) {
  cat("not clean:\n")
  print(bad, row.names = FALSE)
} else {
  cat("not clean: none\n")
}
noresult <- grep("NO RESULT|LOADERROR", lg, value = TRUE)
cat("no result lines:", length(noresult), "\n")
if (nzchar(want_lib)) {
  fl <- list.files(file.path(wt, "dev/dblbar-tests", paste0(tier, "-files")),
                   full.names = TRUE)
  ok <- vapply(fl, function(f) {
    x <- grep("^lib: ", readLines(f, warn = FALSE), value = TRUE)
    length(x) == 1L && grepl(paste0("frmtmb 0.70.0 ", want_lib), x,
                             fixed = TRUE)
  }, NA)
  cat("per-file logs naming the lane core:", sum(ok), "of", length(fl), "\n")
}
