# Reviewer, lane ordinal: per-package totals of dev/ordinal-rev-suite-all,
# from the RESULT lines, plus the lib: line of every log and the failing
# expectations. Output: dev/ordinal-rev-log-suite-sum.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
dir <- file.path(wt, "dev/ordinal-rev-suite-all")
x <- readLines(file.path(dir, "all.log"))
cat(tail(x, 1), "\n")
x <- grep(" RESULT ", x, value = TRUE)
m <- regmatches(x, regexec(paste0("^(\\S+) RESULT \\S+ \\S+ (\\S+) pass=(\\d+) ",
                                  "fail=(\\d+) err=(\\d+) skip=(\\d+) ",
                                  "warn=(\\d+)"), x))
f <- do.call(rbind, lapply(m, function(v) {
  data.frame(pkg = v[2], file = v[3], pass = as.integer(v[4]),
             fail = as.integer(v[5]), err = as.integer(v[6]),
             skip = as.integer(v[7]), warn = as.integer(v[8]))
}))
tot <- aggregate(cbind(files = 1L, pass, fail, err, skip, warn) ~ pkg,
                 data = transform(f, files = 1L), FUN = sum)
print(tot, row.names = FALSE)
cat("\nall:", nrow(f), "files, pass", sum(f$pass), "fail", sum(f$fail),
    "err", sum(f$err), "skip", sum(f$skip), "warn", sum(f$warn), "\n")
cat("\nfiles with a failure, an error, a skip or a warning:\n")
print(f[f$fail + f$err + f$skip + f$warn > 0, ], row.names = FALSE)
logs <- list.files(dir, pattern = "[.]txt$", full.names = TRUE)
libs <- vapply(logs, function(l) {
  z <- grep("^lib:", readLines(l, warn = FALSE), value = TRUE)
  if (length(z)) z[1] else "NO LIB LINE"
}, "")
cat("\nlib: lines (distinct, with counts):\n")
print(table(sub("/[^/ ]+  frmtmb:", "  frmtmb:", libs)))
cat("\nfailure headlines:\n")
for (l in logs) {
  z <- readLines(l, warn = FALSE)
  h <- grep("^── (Failure|Error)", z, value = TRUE)
  if (length(h)) {
    cat(basename(l), "\n")
    for (i in grep("^── (Failure|Error)", z)) cat("   ", z[i], "\n    ", z[i + 1], "\n")
  }
}
