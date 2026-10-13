# Reviewer of lane optrap: dev/optrap-rev-r2.R, warnings and codes per k.
for (b in c("lane", "r7")) {
  x <- read.delim(sprintf("dev/optrap-rev-log/r2-%s.tsv", b))
  cat(b, "
"); print(table(k = x$k, code = x$code, warned = x$nwarn > 0))
}
