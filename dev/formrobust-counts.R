# Generate the test-count block of dev/formrobust-findings.md from the
# per-file logs, so that no count is typed. Usage:
#   Rscript dev/formrobust-counts.R > dev/formrobust-log/counts.md
root <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-log"
read_tier <- function(dir) {
  fs <- list.files(file.path(root, dir), full.names = TRUE)
  rows <- lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    r <- grep("^RESULT", x, value = TRUE)
    lib <- grep("^lib:", x, value = TRUE)
    pkg <- sub("--.*$", "", basename(f))
    if (!length(r)) {
      return(data.frame(pkg = pkg, file = basename(f), pass = NA, fail = NA,
                        err = NA, skip = NA, warn = NA, lib = NA))
    }
    r <- r[length(r)]
    g <- function(k) as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", r))
    data.frame(pkg = pkg, file = strsplit(r, " ")[[1]][2], pass = g("pass"),
               fail = g("fail"), err = g("err"), skip = g("skip"),
               warn = g("warn"),
               lib = if (length(lib)) sub(" frmtmb=.*", "", sub("^lib: ", "",
                                                               lib[1]))
                     else NA)
  })
  do.call(rbind, rows)
}
show_tier <- function(dir, title) {
  d <- read_tier(dir)
  cat("### ", title, " (`dev/formrobust-log/", dir, "`)\n\n", sep = "")
  cat("| package | files | with RESULT | pass | fail | err | skip | warn |\n")
  cat("|---|---|---|---|---|---|---|---|\n")
  for (p in unique(d$pkg)) {
    s <- d[d$pkg == p, ]
    cat(sprintf("| %s | %d | %d | %d | %d | %d | %d | %d |\n", p, nrow(s),
                sum(!is.na(s$pass)), sum(s$pass, na.rm = TRUE),
                sum(s$fail, na.rm = TRUE), sum(s$err, na.rm = TRUE),
                sum(s$skip, na.rm = TRUE), sum(s$warn, na.rm = TRUE)))
  }
  bad <- d[is.na(d$pass) | d$fail > 0 | d$err > 0 | d$warn > 0, ]
  if (nrow(bad)) {
    cat("\nFiles with a failure, an error, a warning or no RESULT line:\n\n")
    for (i in seq_len(nrow(bad))) {
      cat(sprintf("- %s %s: pass=%s fail=%s err=%s skip=%s warn=%s\n",
                  bad$pkg[i], bad$file[i], bad$pass[i], bad$fail[i],
                  bad$err[i], bad$skip[i], bad$warn[i]))
    }
  }
  libs <- unique(d$lib[!is.na(d$lib)])
  cat("\nLibraries the files loaded (package path):\n\n")
  for (l in libs) cat("- ", sum(d$lib == l, na.rm = TRUE), " x ", l, "\n",
                      sep = "")
  cat("\n")
}
args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  parts <- strsplit(a, "=", fixed = TRUE)[[1]]
  show_tier(parts[1], parts[2])
}
