# Generate the suite block of dev/nanse-findings.md from a
# dev/nanse-suite.sh output directory: per package the files, pass,
# fail, err, skip, warn; the lib line of every log; and every firing of
# the standard-error warning by test file.
#   Rscript dev/nanse-suite-sum.R <dir>
d <- commandArgs(TRUE)[1]
txt <- list.files(d, pattern = "[.]txt$", full.names = TRUE)
txt <- txt[basename(txt) != "jobs.txt"]
rows <- lapply(txt, function(f) {
  x <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", x, value = TRUE)
  lib <- grep("^lib: ", x, value = TRUE)
  b <- sub("[.]txt$", "", basename(f))
  num <- function(k) {
    if (!length(r)) return(NA_integer_)
    as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r[1]))
  }
  data.frame(pkg = sub("--.*", "", b), file = sub(".*--", "", b),
             pass = num("pass"), fail = num("fail"), err = num("err"),
             skip = num("skip"), warn = num("warn"),
             lane_core = length(lib) == 1L &&
               grepl("^lib: C:/Users/adf44/source/r/wt-nanse-lib/frmtmb ",
                     lib),
             stringsAsFactors = FALSE)
})
X <- do.call(rbind, rows)
cat("logs:", nrow(X), " with a RESULT line:", sum(!is.na(X$pass)),
    " loading the lane core:", sum(X$lane_core), "\n\n")
for (p in unique(X$pkg)) {
  Y <- X[X$pkg == p, ]
  cat(sprintf("%-16s files %3d pass %5d fail %d err %d skip %d warn %d\n",
              p, nrow(Y), sum(Y$pass, na.rm = TRUE), sum(Y$fail, na.rm = TRUE),
              sum(Y$err, na.rm = TRUE), sum(Y$skip, na.rm = TRUE),
              sum(Y$warn, na.rm = TRUE)))
}
bad <- X[is.na(X$pass) | X$fail > 0 | X$err > 0 | X$warn > 0, ]
if (nrow(bad)) {
  cat("\nnot clean:\n")
  print(bad, row.names = FALSE)
}
sk <- X[!is.na(X$skip) & X$skip > 0, c("pkg", "file", "skip")]
if (nrow(sk)) {
  cat("\nfiles with a skip:\n")
  print(sk, row.names = FALSE)
}
fire <- list.files(d, pattern = "[.]fire$", full.names = TRUE)
n <- vapply(fire, function(f) length(readLines(f, warn = FALSE)), 0L)
cat("\nfirings of the standard-error warning:", sum(n), "in", length(fire),
    "files\n")
for (i in seq_along(fire)) {
  cat(sprintf("  %-45s %d\n", sub("[.]fire$", "", basename(fire[i])), n[i]))
}
