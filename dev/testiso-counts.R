# Summarize one run of dev/release/run-sessions.sh as markdown, from the
# per-session files, for dev/testiso-findings.md.
#
#   Rscript dev/testiso-counts.R <tier>
d <- file.path("dev/release", paste0(commandArgs(TRUE)[1], "-sessions"))
fs <- sort(list.files(d, "[.]txt$", full.names = TRUE))

# The test's own LEAK lines. The runner marks a difference that came with
# a namespace load "[with load of ...]"; the before arm ran an earlier
# runner that did not, so an unmarked line is judged from the LOAD line
# the same file printed, by the same rule: a registry difference is a
# load's when the file loaded an frmtmb namespace, an option or an
# environment variable when it loaded any namespace, and a global or the
# search path is never a load's
own_leaks <- function(x) {
  leak <- grep("^LEAK", x, value = TRUE)
  leak <- leak[!grepl("[[]with load of", leak)]
  load <- grep("^LOAD", x, value = TRUE)
  keep <- vapply(leak, function(l) {
    file <- strsplit(l, " ")[[1]][2]
    ld <- load[startsWith(load, paste0("LOAD ", file, " "))]
    if (!length(ld)) return(TRUE)
    k <- sub("^LEAK [^ ]+ ([^:]+):.*$", "\\1", l)
    if (k %in% c("globals", "search path")) return(TRUE)
    if (k %in% c("options", "envvars")) return(FALSE)
    !grepl(" frmtmb", ld[1], fixed = TRUE)
  }, TRUE)
  leak[keep]
}

cat("From `", d, "/` (", length(fs), " sessions).\n\n", sep = "")
cat("| package | order | files | pass | fail | err | skip | warn |",
    "own LEAK | frmtmb from |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
for (f in fs) {
  x <- readLines(f, warn = FALSE)
  s <- grep("^SESSION", x, value = TRUE)
  lib <- sub("^lib: frmtmb [^ ]+ ", "", grep("^lib: frmtmb ", x,
                                              value = TRUE))
  if (!length(s)) {
    cat("| ", basename(f), " | NO SESSION LINE | | | | | | | | |\n")
    next
  }
  v <- strsplit(s, " ")[[1]]
  kv <- function(k) {
    sub(paste0("^", k, "="), "", grep(paste0("^", k, "="), v, value = TRUE))
  }
  cat("| ", v[2], " | ", v[3], " | ", kv("files"), " | ", kv("pass"),
      " | ", kv("fail"), " | ", kv("err"), " | ", kv("skip"), " | ",
      kv("warn"), " | ", length(own_leaks(x)), " | ", basename(lib),
      " |\n", sep = "")
}
cat("\nFailing and erroring tests (expectation count):\n\n")
none <- TRUE
for (f in fs) {
  x <- readLines(f, warn = FALSE)
  det <- grep("^DETAIL (FAIL|ERROR) ", x, value = TRUE)
  if (!length(det)) next
  none <- FALSE
  tab <- table(sub("^DETAIL (FAIL|ERROR) ([^ ]+) [[]([^]]+)[]].*$",
                   "\\1 \\2 [\\3]", det))
  for (k in names(tab)) {
    cat("- ", sub("[.]txt$", "", basename(f)), ": ", k, " x", tab[[k]],
        "\n", sep = "")
  }
}
if (none) cat("- none\n")
cat("\nLEAK lines a test made itself (not a namespace load):\n\n")
none <- TRUE
for (f in fs) {
  for (l in own_leaks(readLines(f, warn = FALSE))) {
    none <- FALSE
    cat("- ", sub("[.]txt$", "", basename(f)), ": ",
        substr(sub("^LEAK ", "", l), 1, 160), "\n", sep = "")
  }
}
if (none) cat("- none\n")
