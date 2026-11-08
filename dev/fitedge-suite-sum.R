# Sums of the one-file-per-process suite (dev/fitedge-out/fitedge-all.log,
# written by dev/fitedge-run-par.sh), per package, from its RESULT lines.
o <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-out"
x <- readLines(file.path(o, "fitedge-all.log"))
pkg <- NA_character_
rows <- list()
num <- function(l, k) {
  as.integer(sub(paste0(".* ", k, "=([0-9]+).*"), "\\1", paste0(" ", l)))
}
for (l in x) {
  if (startsWith(l, "== ") && endsWith(l, " files ==")) {
    pkg <- strsplit(l, " ", fixed = TRUE)[[1]][2]
  }
  if (startsWith(l, "RESULT ")) {
    rows[[length(rows) + 1L]] <- data.frame(
      pkg = pkg, files = 1L, pass = num(l, "pass"), fail = num(l, "fail"),
      err = num(l, "err"), skip = num(l, "skip"), warn = num(l, "warn"))
  }
}
d <- do.call(rbind, rows)
a <- aggregate(cbind(files, pass, fail, err, skip, warn) ~ pkg, d, sum)
print(a, row.names = FALSE)
cat(sprintf("TOTAL files %d pass %d fail %d err %d skip %d warn %d\n", nrow(d),
            sum(d$pass), sum(d$fail), sum(d$err), sum(d$skip), sum(d$warn)))
cat(tail(x, 1), "\n")
