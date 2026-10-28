# Totals per package of a dev/release/run-par.sh tier log, and the files
# with a failure, an error or a skip.
#   Rscript dev/priordraws-tally.R dev/release/<tier>.log
f <- commandArgs(TRUE)[1]
L <- readLines(f)
pkg <- NA_character_
rows <- list()
for (l in L) {
  if (grepl("^== .* files ==$", l)) pkg <- sub("^== ([^ ]+) .*", "\\1", l)
  if (grepl("^RESULT ", l)) {
    g <- function(k) {
      v <- regmatches(l, regexpr(paste0(k, "=[0-9]+"), l))
      if (length(v)) as.integer(sub(".*=", "", v)) else NA_integer_
    }
    rows[[length(rows) + 1L]] <- data.frame(
      pkg = pkg, file = sub("^RESULT ([^ ]+) .*", "\\1", l),
      pass = g("pass"), fail = g("fail"), err = g("err"), skip = g("skip"),
      warn = g("warn"), loaderror = grepl("LOADERROR", l))
  }
}
d <- do.call(rbind, rows)
agg <- aggregate(cbind(files = 1, pass, fail, err, skip, warn) ~ pkg, d, sum)
print(agg, row.names = FALSE)
cat("TOTAL files", nrow(d), "pass", sum(d$pass), "fail", sum(d$fail),
    "err", sum(d$err), "skip", sum(d$skip), "warn", sum(d$warn),
    "loaderror", sum(d$loaderror), "\n")
bad <- d[d$fail > 0 | d$err > 0 | d$skip > 0 | d$loaderror, ]
if (nrow(bad)) print(bad, row.names = FALSE)
cat(grep(" ran ", L, value = TRUE), "\n")
