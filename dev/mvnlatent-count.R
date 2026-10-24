# Test accounting for the move, generated from the RESULT lines of two
# runs of dev/mvnlatent-par.sh: per package totals, files only in one
# arm, and every file whose counts changed.
#
#   Rscript dev/mvnlatent-count.R <before.log> <after.log>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(f) {
  l <- grep(" RESULT ", readLines(f), value = TRUE)
  get <- function(k) as.numeric(sub(paste0(".* ", k, "=([0-9.]+).*"), "\\1", l))
  data.frame(pkg = sub(" .*", "", l),
             file = sub(".* RESULT ([^ ]+) .*", "\\1", l),
             pass = get("pass"), fail = get("fail"), err = get("err"),
             skip = get("skip"), warn = get("warn"), tests = get("tests"),
             elapsed = get("elapsed"), stringsAsFactors = FALSE)
}
b <- rd(a[1])
x <- rd(a[2])
tot <- function(d) {
  s <- aggregate(cbind(files = 1, pass, fail, err, skip, warn, tests,
                       elapsed) ~ pkg, d, sum)
  s[order(s$pkg), ]
}
cat("BEFORE\n"); print(tot(b), row.names = FALSE)
cat("AFTER\n"); print(tot(x), row.names = FALSE)
kb <- paste(b$pkg, b$file)
kx <- paste(x$pkg, x$file)
cat("\nfiles only before:", paste(setdiff(kb, kx), collapse = ", "), "\n")
cat("files only after:", paste(setdiff(kx, kb), collapse = ", "), "\n")
m <- merge(b, x, by = c("pkg", "file"), all = TRUE,
           suffixes = c(".b", ".a"))
ch <- with(m, is.na(pass.b) | is.na(pass.a) | pass.b != pass.a |
             tests.b != tests.a | fail.b != fail.a | err.b != err.a |
             skip.b != skip.a | warn.b != warn.a)
cat("\nfiles whose counts changed:\n")
print(m[ch, c("pkg", "file", "tests.b", "tests.a", "pass.b", "pass.a",
              "fail.b", "fail.a", "err.b", "err.a", "skip.b", "skip.a",
              "warn.b", "warn.a")], row.names = FALSE)
cat("\nfailing, erroring or warning files after:\n")
print(x[x$fail > 0 | x$err > 0 | x$warn > 0, ], row.names = FALSE)
cat("\nfailing, erroring or warning files before:\n")
print(b[b$fail > 0 | b$err > 0 | b$warn > 0, ], row.names = FALSE)
