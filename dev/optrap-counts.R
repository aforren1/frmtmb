# Lane optrap: the suite table of a dev/ciharden-scan.sh run, generated
# for dev/optrap-findings.md from the run's per-file logs.
#   [OPTRAP_COUNT_LIB=<lib dir name>] Rscript dev/optrap-counts.R <run> [<run> ...]
for (run in commandArgs(TRUE)) {
  dir <- file.path("dev/ciharden-log", paste0(run, "-suite"))
  fs <- list.files(dir, pattern = "[.]txt$", full.names = TRUE)
  jobs <- readLines(file.path("dev/ciharden-log", paste0(run, ".jobs")))
  rows <- lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    r <- grep("^RESULT", x, value = TRUE)
    l <- grep("^lib: ", x, value = TRUE)
    num <- function(k) {
      m <- regmatches(r, regexpr(paste0(k, "=[0-9]+"), r))
      if (length(m)) as.integer(sub(".*=", "", m[length(m)])) else NA
    }
    data.frame(file = basename(f),
               pkg = sub("--.*", "", basename(f)),
               result = length(r) > 0,
               lane_lib = length(l) > 0 &&
                 grepl(Sys.getenv("OPTRAP_COUNT_LIB", "wt-optrap-lib"), l[1], fixed = TRUE),
               pass = num("pass"), fail = num("fail"), err = num("err"),
               skip = num("skip"), warn = num("warn"),
               stringsAsFactors = FALSE)
  })
  X <- do.call(rbind, rows)
  cat(sprintf("\n### %s\n\n", run))
  cat(sprintf(paste0("%d jobs; %d files with a RESULT line; %d whose ",
                     "`lib:` line names the lane's library.\n\n"),
              length(jobs), sum(X$result), sum(X$lane_lib)))
  cat("| package | files | pass | fail | error | skip | warn |\n")
  cat("|---|---|---|---|---|---|---|\n")
  for (p in unique(X$pkg)) {
    x <- X[X$pkg == p, ]
    cat(sprintf("| %s | %d | %d | %d | %d | %d | %d |\n", p, nrow(x),
                sum(x$pass, na.rm = TRUE), sum(x$fail, na.rm = TRUE),
                sum(x$err, na.rm = TRUE), sum(x$skip, na.rm = TRUE),
                sum(x$warn, na.rm = TRUE)))
  }
  cat(sprintf("| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n",
              nrow(X), sum(X$pass, na.rm = TRUE), sum(X$fail, na.rm = TRUE),
              sum(X$err, na.rm = TRUE), sum(X$skip, na.rm = TRUE),
              sum(X$warn, na.rm = TRUE)))
  bad <- X[!X$result | X$fail > 0 | X$err > 0 | X$warn > 0, ]
  cat(sprintf("\nFiles with a failure, an error, an escaped warning or no result: %d\n\n",
              nrow(bad)))
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("- %s: pass %s fail %s err %s warn %s\n", bad$file[i],
                bad$pass[i], bad$fail[i], bad$err[i], bad$warn[i]))
    d <- grep("^DETAIL\t(FAIL|ERROR|WARN)",
              readLines(file.path(dir, bad$file[i]), warn = FALSE),
              value = TRUE)
    for (z in d) cat("    -", substr(gsub("\t", " ", z), 8, 260), "\n")
  }
}
