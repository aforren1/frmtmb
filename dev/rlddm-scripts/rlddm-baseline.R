# Join a suite run against dev/suite-baseline.tsv, per file.
#
#   Rscript dev/rlddm-scripts/rlddm-baseline.R <pkg> <suite-log>
#
# The baseline file is what it exists for. A file this lane did not
# touch that moves off its baseline count is a regression; a file it did
# touch is expected to move and by how much is stated in
# dev/rlddm-findings.md.
#
# A file whose runner produced no RESULT line is a HOLE, not a zero, and
# is printed as one. A runner that summed it as zero would report a
# clean total for a file that never ran, which is one of the five broken
# measurements dev/round-handoff.md records.

args <- commandArgs(trailingOnly = TRUE)
pkg <- args[[1L]]
log <- args[[2L]]

base <- utils::read.delim(
  "dev/suite-baseline.tsv", header = FALSE,
  col.names = c("pkg", "file", "pass", "skip"),
  stringsAsFactors = FALSE)
base <- base[base$pkg == pkg, , drop = FALSE]

ln <- grep("^RESULT ", readLines(log), value = TRUE)
grab <- function(x, k) {
  v <- rep(NA_real_, length(x))
  hit <- grepl(paste0("\\b", k, "=[0-9-]+"), x)
  v[hit] <- as.numeric(sub(paste0(".*\\b", k, "=(-?[0-9]+).*"), "\\1",
                           x[hit]))
  v
}
got <- data.frame(
  file = sub(".*file=([^ ]+).*", "\\1", ln),
  pass = grab(ln, "pass"), fail = grab(ln, "fail"),
  err = grab(ln, "err"), skip = grab(ln, "skip"),
  exit = grab(ln, "exit"), stringsAsFactors = FALSE)
got$hole <- is.na(got$pass)

m <- merge(got, base[, c("file", "pass", "skip")], by = "file",
           all.x = TRUE, suffixes = c("", "_base"))
m$delta <- m$pass - m$pass_base
n <- function(z) if (is.na(z)) "-" else format(z)
cat(sprintf("%-28s %6s %6s %5s %5s %6s %6s %s\n", "file", "pass",
            "base", "fail", "err", "skip", "exit", "note"))
for (i in order(m$file)) {
  note <- if (isTRUE(m$hole[i])) {
    "NO RESULT LINE"
  } else if (is.na(m$pass_base[i])) {
    "NEW"
  } else if (m$delta[i] == 0) {
    "at baseline"
  } else {
    sprintf("%+d", m$delta[i])
  }
  cat(sprintf("%-28s %6s %6s %5s %5s %6s %6s %s\n", m$file[i],
              n(m$pass[i]),
              if (is.na(m$pass_base[i])) "-" else format(m$pass_base[i]),
              n(m$fail[i]), n(m$err[i]), n(m$skip[i]), n(m$exit[i]),
              note))
}
cat("\nfiles                    ", nrow(got), "\n")
cat("files with NO RESULT line", sum(got$hole), "\n")
cat("assertions               ", sum(got$pass, na.rm = TRUE), "\n")
cat("failures                 ", sum(got$fail, na.rm = TRUE), "\n")
cat("errors                   ", sum(got$err, na.rm = TRUE), "\n")
cat("skips                    ", sum(got$skip, na.rm = TRUE), "\n")
cat("non-zero exits           ", sum(got$exit != 0, na.rm = TRUE), "\n")
cat("off baseline, not NEW    ",
    sum(!is.na(m$pass_base) & !is.na(m$delta) & m$delta != 0), "\n")
