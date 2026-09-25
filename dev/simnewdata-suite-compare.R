# Sum a tier log and compare it file by file with dev/suite-baseline.tsv.
#   Rscript dev/simnewdata-suite-compare.R dev/simnewdata-log/suite.log
f <- commandArgs(trailingOnly = TRUE)[1]
x <- readLines(f)
pkg <- NA_character_
rows <- list()
for (l in x) {
  if (grepl("^== .* files ==$", l)) pkg <- sub("^== (\\S+) .*$", "\\1", l)
  if (!startsWith(l, "RESULT ")) next
  m <- regmatches(l, regexec(
    "^RESULT (\\S+) pass=([0-9]+) fail=([0-9]+) err=([0-9]+) skip=([0-9]+)", l))[[1]]
  if (!length(m)) {
    rows[[length(rows) + 1L]] <- data.frame(pkg = pkg, file = l, pass = NA,
                                            fail = NA, err = NA, skip = NA)
    next
  }
  rows[[length(rows) + 1L]] <- data.frame(pkg = pkg, file = m[2],
    pass = as.integer(m[3]), fail = as.integer(m[4]), err = as.integer(m[5]),
    skip = as.integer(m[6]))
}
r <- do.call(rbind, rows)
cat("files", nrow(r), "pass", sum(r$pass, na.rm = TRUE), "fail",
    sum(r$fail, na.rm = TRUE), "err", sum(r$err, na.rm = TRUE), "skip",
    sum(r$skip, na.rm = TRUE), "unparsed", sum(is.na(r$pass)), "\n")
bad <- r[is.na(r$pass) | r$fail > 0 | r$err > 0, ]
if (nrow(bad)) print(bad)
b <- utils::read.delim("dev/suite-baseline.tsv", header = FALSE,
                       col.names = c("pkg", "file", "pass", "skip"))
m <- merge(r, b, by = c("pkg", "file"), all = TRUE,
           suffixes = c("", ".base"))
cat("\nfiles whose count moved against the 0.62.0 baseline:\n")
mv <- m[is.na(m$pass) | is.na(m$pass.base) | m$pass != m$pass.base |
          m$skip != m$skip.base, ]
print(mv[, c("pkg", "file", "pass", "pass.base", "skip", "skip.base")],
      row.names = FALSE)
