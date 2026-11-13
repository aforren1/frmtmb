# Lane nlstop: per-package totals of a run-par.sh tier from its per-file
# logs, with the files whose lib: line does not name the lane's core, and
# every file with a failure, error, escaped warning or no RESULT line.
# Usage: Rscript dev/nlstop-suitesum.R <dir of per-file logs> [lib substring]
a <- commandArgs(TRUE)
dir <- a[1]
want <- if (length(a) > 1) a[2] else "wt-nlstop-lib"
fs <- list.files(dir, "[.]txt$", full.names = TRUE)
rows <- lapply(fs, function(f) {
  l <- readLines(f, warn = FALSE)
  r <- tail(grep("^RESULT", l, value = TRUE), 1)
  lib <- grep("^lib:", l, value = TRUE)[1]
  pkg <- sub("--.*$", "", basename(f))
  num <- function(k) {
    if (!length(r)) return(NA_integer_)
    as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
  }
  data.frame(pkg = pkg, file = basename(f), pass = num("pass"),
             fail = num("fail"), err = num("err"), skip = num("skip"),
             warn = num("warn"),
             core_ok = !is.na(lib) && grepl(paste0("frmtmb 0.81.0 [^;]*", want),
                                            sub("^.*; ", "", lib)))
})
X <- do.call(rbind, rows)
cat("files", nrow(X), "with RESULT", sum(!is.na(X$pass)),
    "with the lane's core on the lib: line", sum(X$core_ok), "\n")
agg <- aggregate(cbind(files = 1, pass, fail, err, skip, warn) ~ pkg,
                 data = transform(X, files = 1), FUN = sum, na.action = na.pass)
print(agg, row.names = FALSE)
cat("total pass", sum(X$pass, na.rm = TRUE), "fail", sum(X$fail, na.rm = TRUE),
    "err", sum(X$err, na.rm = TRUE), "skip", sum(X$skip, na.rm = TRUE),
    "warn", sum(X$warn, na.rm = TRUE), "\n")
bad <- X[is.na(X$pass) | X$fail > 0 | X$err > 0 | X$warn > 0 | !X$core_ok, ]
cat("files with a failure, error, warning, no RESULT or another core:",
    nrow(bad), "\n")
if (nrow(bad)) print(bad, row.names = FALSE)
sk <- X[!is.na(X$skip) & X$skip > 0, c("pkg", "file", "skip")]
cat("files with a skip:", nrow(sk), "\n")
writeLines(paste(sk$pkg, sk$file), file.path(dirname(dir),
                                             paste0(basename(dir), "-skips.lst")))
