# Lane postfix: the rows of a brms-suite record (brmsport-bin2-run.R or
# the bin-1 recorder) whose outcome differs between two records.
#   Rscript dev/postfix-bin2-diff.R <before.tsv> <after.tsv>
a <- commandArgs(TRUE)
rd <- function(f) {
  d <- utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE,
                         quote = "")
  d <- d[d$V1 %in% c("assert", "setup"), ]
  d$key <- paste(d$V1, d$V3)
  d
}
b <- rd(a[1]); w <- rd(a[2])
m <- merge(b[c("key", "V4", "V5", "V7")], w[c("key", "V5", "V7")],
           by = "key", suffixes = c(".b", ".a"), all = TRUE)
ch <- m[!(m$V5.b == m$V5.a) %in% TRUE, ]
cat("rows:", nrow(m), "; outcome changed:", nrow(ch), "\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf(paste0("%-28s verdict %-10s held/ok %s -> %s\n",
                     "    before: %s\n    after:  %s\n"),
              ch$key[i], ch$V4[i], ch$V5.b[i], ch$V5.a[i],
              substr(ch$V7.b[i], 1, 160), substr(ch$V7.a[i], 1, 160)))
}
