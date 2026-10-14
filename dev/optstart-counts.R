# Lane optstart: per-package counts of a scan (dev/ciharden-scan.sh,
# dev/ciharden-log/<run>.sum and <run>-suite/): files, pass, fail,
# error, skip, warn, how many files' lib: line names the expected
# library, and every file with a failure, an error, an escaped warning
# or no result line.
#   Rscript dev/optstart-counts.R <run> <expected frmtmb library>
a <- commandArgs(TRUE)
fs <- Sys.glob(file.path("dev/ciharden-log", paste0(a[1], "-suite"), "*.txt"))
rows <- lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", x, value = TRUE)
  lib <- grep("^lib: ", x, value = TRUE)
  pkg <- sub("--.*", "", basename(f))
  num <- function(k) {
    if (!length(r)) return(NA_integer_)
    as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r[length(r)]))
  }
  data.frame(pkg = pkg, file = sub(".*--", "", basename(f)),
             result = length(r) > 0,
             libok = length(lib) > 0 && grepl(a[2], lib[1], fixed = TRUE),
             pass = num("pass"), fail = num("fail"), err = num("err"),
             skip = num("skip"), warn = num("warn"))
})
X <- do.call(rbind, rows)
cat(sprintf("%d files; %d with a RESULT line; %d whose lib: line names %s\n\n",
            nrow(X), sum(X$result), sum(X$libok), a[2]))
cat("| package | files | pass | fail | error | skip | warn |\n")
cat("|---|---|---|---|---|---|---|\n")
for (p in sort(unique(X$pkg))) {
  x <- X[X$pkg == p, ]
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", p, nrow(x),
              sum(x$pass, na.rm = TRUE), sum(x$fail, na.rm = TRUE),
              sum(x$err, na.rm = TRUE), sum(x$skip, na.rm = TRUE),
              sum(x$warn, na.rm = TRUE)))
}
cat(sprintf("| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n",
            nrow(X), sum(X$pass, na.rm = TRUE), sum(X$fail, na.rm = TRUE),
            sum(X$err, na.rm = TRUE), sum(X$skip, na.rm = TRUE),
            sum(X$warn, na.rm = TRUE)))
bad <- X[!X$result | X$fail > 0 | X$err > 0 | X$warn > 0, ]
cat("\nfiles with a failure, an error, an escaped warning or no result:",
    nrow(bad), "\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf("- %s/%s: pass %s fail %s err %s skip %s warn %s\n",
              bad$pkg[i], bad$file[i], bad$pass[i], bad$fail[i], bad$err[i],
              bad$skip[i], bad$warn[i]))
}
