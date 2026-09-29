# Lane wt-resmooth. Generate the suite counts block from the per-file
# result lines, so no count in dev/resmooth-findings.md is typed.
#   Rscript dev/resmooth-summarise.R > dev/resmooth-suite-counts.txt
read_dir <- function(dir, label) {
  fs <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  rows <- list()
  for (f in fs) {
    x <- readLines(f, warn = FALSE)
    ln <- grep("^RESULT ", x, value = TRUE)
    if (!length(ln)) {
      rows[[length(rows) + 1L]] <- data.frame(
        file = basename(f), pass = NA_integer_, fail = NA_integer_,
        error = NA_integer_, skip = NA_integer_, warn = NA_integer_,
        note = "NO RESULT LINE")
      next
    }
    g <- function(k) {
      as.integer(sub(paste0(".*\\b", k, "=([0-9]+).*"), "\\1", ln[[1L]]))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      file = basename(f), pass = g("pass"), fail = g("fail"),
      error = g("error"), skip = g("skip"), warn = g("warn"), note = "")
  }
  d <- do.call(rbind, rows)
  # wrapped so the pasted block stays inside 80 columns, which the
  # house style asks of markdown too
  cat(sprintf("%s: %d files\n  pass=%d fail=%d error=%d skip=%d warn=%d\n",
              label, nrow(d), sum(d$pass, na.rm = TRUE),
              sum(d$fail, na.rm = TRUE), sum(d$error, na.rm = TRUE),
              sum(d$skip, na.rm = TRUE), sum(d$warn, na.rm = TRUE)))
  bad <- d[is.na(d$pass) | d$fail > 0 | d$error > 0, , drop = FALSE]
  if (nrow(bad)) {
    cat("  files with a failure, an error or no result line:\n")
    for (i in seq_len(nrow(bad))) {
      cat(sprintf("  - %s pass=%s fail=%s error=%s skip=%s %s\n",
                  bad$file[i], bad$pass[i], bad$fail[i], bad$error[i],
                  bad$skip[i], bad$note[i]))
    }
  } else {
    cat("  no file failed, errored or failed to report\n")
  }
  sk <- d[!is.na(d$skip) & d$skip > 0, , drop = FALSE]
  cat(sprintf("  files with a skip: %d\n", nrow(sk)))
  for (i in seq_len(nrow(sk))) {
    cat(sprintf("    %s skip=%d\n", sk$file[i], sk$skip[i]))
  }
  invisible(d)
}
# A count from a run that has not finished is not a count: a redirect
# creates the output file before the R process writes anything, so the
# file count reaches 179 long before the suite does. The driver's own
# last line is the positive condition, and it was OBSERVED to be absent
# while four brms files were still running.
logs <- strsplit(Sys.getenv(
  "RESMOOTH_DRIVERLOGS",
  "dev/resmooth-suite.driverlog,dev/resmooth-ext.driverlog"), ",")[[1L]]
for (lg in logs) {
  txt <- if (file.exists(lg)) readLines(lg, warn = FALSE) else character(0)
  if (!any(grepl("-DONE$", txt))) {
    cat("REFUSING TO COUNT:", lg, "does not say DONE yet\n")
    quit(status = 1L)
  }
}
cat("== GENERATED BLOCK, paste verbatim ==\n")
read_dir(Sys.getenv("RESMOOTH_SUITE_DIR", "dev/resmooth-suite"),
         "core frmtmb suite")
read_dir(Sys.getenv("RESMOOTH_EXT_DIR", "dev/resmooth-ext"),
         "frmtmb.sample + frmtmb.spline suites")
cat("== end ==\n")
