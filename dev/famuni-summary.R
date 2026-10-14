# Summarize a famuni suite directory from its per-file logs: the totals,
# the libraries loaded, the files with a failure, an error or no RESULT
# line, and one line per file.
# Usage: Rscript dev/famuni-summary.R <dir>   (writes <dir>/SUMMARY.txt)
d <- commandArgs(TRUE)[1]
fs <- list.files(d, pattern = "--.*[.]txt$", full.names = TRUE)
jobs <- length(readLines(file.path(d, "JOBS.txt")))
keys <- c("tests", "pass", "fail", "error", "skip", "warn")
tab <- t(vapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  r <- grep("^RESULT", x, value = TRUE)
  if (!length(r)) return(rep(NA_real_, length(keys)))
  vapply(keys, function(k) {
    as.numeric(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r[1]))
  }, 0)
}, numeric(length(keys))))
colnames(tab) <- keys
libs <- vapply(fs, function(f) {
  l <- grep("^lib:", readLines(f, warn = FALSE), value = TRUE)
  if (length(l)) sub("^lib: (\\S+).*$", "\\1", l[1]) else NA_character_
}, "")
nm <- sub("[.]txt$", "", basename(fs))
bad <- is.na(tab[, "fail"]) | tab[, "fail"] > 0 | tab[, "error"] > 0
lt <- table(libs, useNA = "ifany")
out <- c(sprintf("jobs %d, logs %d, with a RESULT line %d", jobs, length(fs),
                 sum(!is.na(tab[, "tests"]))),
         paste("libraries loaded:",
               paste(names(lt), lt, sep = " x", collapse = "; ")),
         paste("totals:", paste(keys, colSums(tab, na.rm = TRUE),
                                collapse = " ")),
         "files with a failure, an error or no RESULT:",
         if (any(bad)) paste(" ", nm[bad], apply(tab[bad, , drop = FALSE], 1,
                                               paste, collapse = " ")) else
           "  none",
         "per file (tests pass fail error skip warn):",
         paste(" ", nm, apply(tab, 1, paste, collapse = " ")))
writeLines(out, file.path(d, "SUMMARY.txt"))
cat(out[1:(5 + sum(bad))], sep = "\n")
