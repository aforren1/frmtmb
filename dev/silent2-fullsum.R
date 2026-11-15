# Lane silent2: per-package totals of a run-par.sh tier, generated from
# its per-file logs, with the lib: lines checked against the lane library.
#   Rscript dev/silent2-fullsum.R <tier>
tier <- commandArgs(TRUE)[1]
dir <- file.path("dev/release", paste0(tier, "-files"))
fs <- list.files(dir, pattern = "[.]txt$", full.names = TRUE)
pkg <- sub("--.*$", "", basename(fs))
num <- function(line, key) {
  as.integer(sub(paste0(".*", key, "=([0-9]+).*"), "\\1", line))
}
rows <- lapply(seq_along(fs), function(i) {
  x <- readLines(fs[i], warn = FALSE)
  r <- utils::tail(grep("^RESULT", x, value = TRUE), 1L)
  lib <- grep("^lib:|frmtmb from|library", x, value = TRUE)
  data.frame(pkg = pkg[i], file = basename(fs[i]), has = length(r) == 1L,
             pass = if (length(r)) num(r, "pass") else NA,
             fail = if (length(r)) num(r, "fail") else NA,
             err = if (length(r)) num(r, "err") else NA,
             skip = if (length(r)) num(r, "skip") else NA,
             warn = if (length(r)) num(r, "warn") else NA,
             lanelib = any(grepl("wt-silent2-lib", x, fixed = TRUE)))
})
R <- do.call(rbind, rows)
for (p in unique(R$pkg)) {
  s <- R[R$pkg == p, ]
  cat(sprintf("%s files=%d pass=%d fail=%d err=%d skip=%d warn=%d lanelib=%d\n",
              p, nrow(s), sum(s$pass, na.rm = TRUE), sum(s$fail, na.rm = TRUE),
              sum(s$err, na.rm = TRUE), sum(s$skip, na.rm = TRUE),
              sum(s$warn, na.rm = TRUE), sum(s$lanelib)))
}
cat("no RESULT line:", sum(!R$has), "\n")
bad <- R[!R$has | (R$fail + R$err + R$warn) > 0, ]
if (nrow(bad)) print(bad, row.names = FALSE)
sk <- R[(R$skip > 0) %in% TRUE, c("pkg", "file", "skip")]
cat("files with skips:\n")
print(sk, row.names = FALSE)
