# Per-package totals of a tier's per-file logs, and the check that each
# log's lib: line names the lane library for frmtmb.
# Usage: Rscript dev/sigmainn-sum.R <dir of per-file logs>
d <- commandArgs(TRUE)[1]
fs <- list.files(d, pattern = "[.]txt$", full.names = TRUE)
rows <- lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  pkg <- sub("--.*$", "", basename(f))
  r <- grep("^RESULT", x, value = TRUE)
  lib <- grep("^lib:", x, value = TRUE)
  n <- function(k) {
    if (!length(r)) return(NA_integer_)
    as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", tail(r, 1)))
  }
  data.frame(pkg = pkg, file = basename(f), has = length(r) > 0 &&
               !grepl("LOADERROR", tail(r, 1)),
             pass = n("pass"), fail = n("fail"), err = n("err"),
             skip = n("skip"), warn = n("warn"),
             mylib = length(lib) > 0 &&
               grepl("frmtmb 0.76.0 C:/Users/adf44/source/r/wt-sigmainn-lib",
                     lib[1], fixed = TRUE))
})
t <- do.call(rbind, rows)
for (p in unique(t$pkg)) {
  s <- t[t$pkg == p, ]
  cat(sprintf("%-16s files=%d pass=%d fail=%d err=%d skip=%d warn=%d\n", p,
              nrow(s), sum(s$pass, na.rm = TRUE), sum(s$fail, na.rm = TRUE),
              sum(s$err, na.rm = TRUE), sum(s$skip, na.rm = TRUE),
              sum(s$warn, na.rm = TRUE)))
}
cat("logs without a RESULT line:", sum(!t$has), "\n")
cat("logs whose frmtmb is not the lane lib:", sum(!t$mylib), "\n")
bad <- t[!t$has | (t$fail %in% NA) | t$fail > 0 | t$err > 0 | t$warn > 0, ]
if (nrow(bad)) print(bad[, c("pkg", "file", "pass", "fail", "err", "warn")])
