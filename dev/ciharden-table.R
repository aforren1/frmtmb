# The before/after table of dev/ciharden-findings.md: for every test
# file the scan found configuration-dependent on the base (b-*) runs,
# its RESULT counts under each configuration, before (base tree on
# rellib-r6) and after (this worktree on the lane library), plus the
# suite totals of every run. Printed as markdown; paste it verbatim.
#   Rscript dev/ciharden-table.R
logd <- "dev/ciharden-log"
cfg <- c("ref", "ob0.3.26-blas-t4", "ob0.3.26-lapack-t4",
         "ob0.3.32-lapack-t4", "ob0.3.32-lapack-t1")
res_of <- function(run, key) {
  f <- file.path(logd, paste0(run, "-suite"), paste0(key, ".txt"))
  if (!file.exists(f)) return("-")
  l <- grep("^RESULT ", readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return("no result")
  m <- regmatches(l, regexpr("pass=.*$", l))
  v <- as.integer(sub(".*=", "", strsplit(m, " ")[[1]]))
  paste(v, collapse = "/")
}
d <- utils::read.delim(file.path(logd, "b-diff.tsv"), check.names = FALSE)
keys <- unique(d$file)
cat("Counts are pass/fail/error/skip/warn.\n\n")
cat("| File | arm |", paste(cfg, collapse = " | "), "|\n")
cat("|---|---|", paste(rep("---", length(cfg)), collapse = "|"), "|\n")
for (k in keys) {
  for (arm in c("b", "a")) {
    v <- vapply(cfg, function(c) res_of(paste0(arm, "-", c), k), "")
    cat("|", if (arm == "b") sub("--", " ", k) else "", "|",
        if (arm == "b") "before" else "after", "|",
        paste(v, collapse = " | "), "|\n")
  }
}
cat("\nSuite totals, every file of the eight suites:\n\n")
cat("| Run | files with a RESULT | pass | fail | error | skip | warn |\n")
cat("|---|---|---|---|---|---|---|\n")
runs <- c(paste0("b-", c(cfg, "ref2")), paste0("a-", cfg))
for (r in runs) {
  s <- file.path(logd, paste0(r, ".sum"))
  dd <- file.path(logd, paste0(r, "-suite"))
  if (!dir.exists(dd)) next
  fs <- list.files(dd, "[.]txt$", full.names = TRUE)
  tot <- integer(5)
  n <- 0L
  for (f in fs) {
    l <- grep("^RESULT ", readLines(f, warn = FALSE), value = TRUE)
    if (!length(l)) next
    m <- regmatches(l, regexpr("pass=.*$", l))
    if (!length(m)) next
    n <- n + 1L
    tot <- tot + as.integer(sub(".*=", "", strsplit(m, " ")[[1]]))
  }
  cat(sprintf("| %s | %d of %d | %s |\n", r, n, length(fs),
              paste(tot, collapse = " | ")))
}
