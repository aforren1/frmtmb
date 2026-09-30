# Reviewer: counts per package from the per-file logs of a tier.
# Usage: Rscript formrobust-rev-counts.R <tier>
tier <- commandArgs(trailingOnly = TRUE)[1]
dir <- file.path("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-log",
                 paste0(tier, "-files"))
fs <- list.files(dir, full.names = TRUE)
rows <- lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  pkg <- sub("--.*$", "", basename(f))
  r <- grep("^RESULT", x, value = TRUE)
  lib <- sub("^lib: (\\S+) .*$", "\\1", grep("^lib: ", x, value = TRUE))
  num <- function(k) {
    if (!length(r)) return(NA_integer_)
    as.integer(sub(paste0("^.* ", k, "=([0-9]+).*$"), "\\1", tail(r, 1)))
  }
  data.frame(pkg = pkg, file = basename(f), result = length(r) > 0,
             pass = num("pass"), fail = num("fail"), err = num("err"),
             skip = num("skip"), warn = num("warn"),
             lib = if (length(lib)) lib[1] else NA_character_)
})
d <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(d, d$pkg), function(s) {
  data.frame(pkg = s$pkg[1], files = nrow(s), with_result = sum(s$result),
             pass = sum(s$pass, na.rm = TRUE), fail = sum(s$fail, na.rm = TRUE),
             err = sum(s$err, na.rm = TRUE), skip = sum(s$skip, na.rm = TRUE),
             warn = sum(s$warn, na.rm = TRUE),
             lib = paste(unique(dirname(s$lib)), collapse = ","))
}))
cat("| package | files | with RESULT | pass | fail | err | skip | warn | library |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(agg))) {
  with(agg[i, ], cat(sprintf("| %s | %d | %d | %d | %d | %d | %d | %d | %s |\n",
                             pkg, files, with_result, pass, fail, err, skip,
                             warn, lib)))
}
cat(sprintf("| total | %d | %d | %d | %d | %d | %d | %d | |\n", nrow(d),
            sum(d$result), sum(d$pass, na.rm = TRUE),
            sum(d$fail, na.rm = TRUE), sum(d$err, na.rm = TRUE),
            sum(d$skip, na.rm = TRUE), sum(d$warn, na.rm = TRUE)))
bad <- d[!d$result | d$fail > 0 | d$err > 0 | d$warn > 0, ]
if (nrow(bad)) {
  cat("\nFiles with a failure, error, warning or no RESULT:\n")
  for (i in seq_len(nrow(bad))) with(bad[i, ], cat(sprintf(
    "- %s %s: pass=%s fail=%s err=%s skip=%s warn=%s\n", pkg, file, pass,
    fail, err, skip, warn)))
}
sk <- d[!is.na(d$skip) & d$skip > 0, ]
cat("\nFiles with skips:", nrow(sk), "\n")
