# The counts of a run-par.sh tier, per package, from its per-file logs: files,
# pass, fail, error, skip, warn, the library each log says it read, and every
# file with a failure, error, skip or escaped warning.
#   Rscript dev/predread-suite-table.R <tier>
tier <- commandArgs(TRUE)[1]
d <- file.path("dev/release", paste0(tier, "-files"))
fs <- list.files(d, full.names = TRUE)
rows <- lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", x, value = TRUE)
  lib <- sub(".*; frmtmb [^ ]* ", "", grep("^lib:", x, value = TRUE)[1L])
  pkg <- sub("--.*", "", basename(f))
  if (!length(r)) {
    return(data.frame(pkg = pkg, file = basename(f), pass = NA, fail = NA,
                      err = NA, skip = NA, warn = NA, lib = lib))
  }
  g <- function(k) as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r[1L]))
  data.frame(pkg = pkg, file = sub("RESULT ([^ ]+).*", "\\1", r[1L]),
             pass = g("pass"), fail = g("fail"), err = g("err"),
             skip = g("skip"), warn = g("warn"), lib = lib)
})
a <- do.call(rbind, rows)
cat("| package | files | pass | fail | error | skip | warn | core read from |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (p in unique(a$pkg)) {
  s <- a[a$pkg == p, ]
  cat(sprintf("| %s | %d | %d | %d | %d | %d | %d | %s |\n", p, nrow(s),
              sum(s$pass, na.rm = TRUE), sum(s$fail, na.rm = TRUE),
              sum(s$err, na.rm = TRUE), sum(s$skip, na.rm = TRUE),
              sum(s$warn, na.rm = TRUE), paste(unique(s$lib), collapse = ", ")))
}
cat("\nfiles with no RESULT:", sum(is.na(a$pass)), "\n\n")
bad <- a[is.na(a$pass) | a$fail > 0 | a$err > 0 | a$warn > 0, ]
for (i in seq_len(nrow(bad))) {
  cat(sprintf("- %s %s: pass %s fail %s err %s skip %s warn %s\n", bad$pkg[i],
              bad$file[i], bad$pass[i], bad$fail[i], bad$err[i], bad$skip[i],
              bad$warn[i]))
}
sk <- a[a$skip %in% seq_len(1e6), ]
cat("\nfiles with a skip:", nrow(sk), "\n")
for (i in seq_len(nrow(sk))) {
  cat(sprintf("- %s %s: skip %d\n", sk$pkg[i], sk$file[i], sk$skip[i]))
}
