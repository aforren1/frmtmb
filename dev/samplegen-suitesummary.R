# Counts for a directory of dev/samplegen-runtests*.R outputs, summed
# from the files; a file without a FILE line is reported as no result.
#   Rscript dev/samplegen-suitesummary.R <outdir> <label>
av <- commandArgs(trailingOnly = TRUE)
fs <- list.files(av[1], pattern = "[.]txt$", full.names = TRUE)
num <- function(ln, key) {
  x <- grep(paste0("^", key, " "), ln, value = TRUE)
  if (!length(x)) NA_real_ else as.numeric(trimws(sub(key, "", x[1])))
}
tab <- t(vapply(fs, function(f) {
  ln <- readLines(f, warn = FALSE)
  c(result = any(grepl("^FILE ", ln)), blocks = num(ln, "BLOCKS"),
    assert = num(ln, "ASSERT"), pass = num(ln, "PASS"),
    fail = num(ln, "FAIL"), error = num(ln, "ERROR"),
    skip = num(ln, "SKIP"), warn = num(ln, "WARN"))
}, numeric(8)))
cat("```\n== suite:", av[2], "==\n")
cat(sprintf("%-26s %d of %d\n", "test files with a result",
            sum(tab[, "result"]), nrow(tab)))
for (k in c("blocks", "assert", "pass", "fail", "error", "skip", "warn")) {
  cat(sprintf("%-26s %d\n", k, as.integer(sum(tab[, k], na.rm = TRUE))))
}
cat("\nper file: blocks / assertions / fail / error / skip\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  %-40s %4d %5d %3d %3d %3d\n",
              sub("[.]txt$", "", basename(rownames(tab)[i])),
              as.integer(tab[i, "blocks"]), as.integer(tab[i, "assert"]),
              as.integer(tab[i, "fail"]), as.integer(tab[i, "error"]),
              as.integer(tab[i, "skip"])))
}
bad <- rownames(tab)[!tab[, "result"] | tab[, "fail"] > 0 |
                       tab[, "error"] > 0]
cat(if (length(bad)) paste("\nNOT CLEAN:", basename(bad)) else
  "\nno file reported a failure or an error", "\n```\n")
