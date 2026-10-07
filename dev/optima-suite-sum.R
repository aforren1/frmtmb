# Lane optima: per-package totals of a run-par tier, from its per-file
# logs: files, pass, fail, error, skip, escaped warnings, and every file
# that is not clean, with its first DETAIL lines.
#   Rscript dev/optima-suite-sum.R <tier>
tier <- commandArgs(TRUE)[1]
fs <- Sys.glob(file.path("dev/release", paste0(tier, "-files"), "*.txt"))
rows <- lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  r <- utils::tail(grep("^RESULT", x, value = TRUE), 1)
  pkg <- sub("--.*", "", basename(f))
  if (!length(r)) {
    return(data.frame(pkg = pkg, file = basename(f), pass = NA, fail = NA,
                      err = NA, skip = NA, warn = NA, ok = FALSE))
  }
  g <- function(k) as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r))
  data.frame(pkg = pkg, file = basename(f), pass = g("pass"),
             fail = g("fail"), err = g("err"), skip = g("skip"),
             warn = g("warn"), ok = TRUE)
})
X <- do.call(rbind, rows)
cat("tier", tier, ":", nrow(X), "files,", sum(!X$ok), "without a RESULT line\n")
agg <- stats::aggregate(cbind(files = 1, pass, fail, err, skip, warn) ~ pkg,
                        data = X[X$ok, ], FUN = sum)
print(agg, row.names = FALSE)
bad <- X[!X$ok | X$fail > 0 | X$err > 0 | X$warn > 0 | X$skip > 0, ]
for (i in seq_len(nrow(bad))) {
  b <- bad[i, ]
  cat(sprintf("-- %s %s: pass %s fail %s err %s skip %s warn %s\n", b$pkg,
              b$file, b$pass, b$fail, b$err, b$skip, b$warn))
  x <- readLines(file.path("dev/release", paste0(tier, "-files"), b$file),
                 warn = FALSE)
  d <- grep("^DETAIL (FAIL|ERROR|WARN)", x, value = TRUE)
  if (length(d)) cat(paste0("   ", substr(utils::head(d, 2), 1, 160)),
                     sep = "\n")
}
