# The bin-2 rows whose outcome moved between the 0.72.0 build and this
# lane's, on the same regenerated distributions and log_lik blocks:
# held (column 5) by kind, id and the id's occurrence
rd <- function(f) {
  x <- read.delim(f, header = FALSE, quote = "", colClasses = "character")
  x$key <- paste(x$V1, x$V3, ave(seq_len(nrow(x)), x$V1, x$V3, FUN = seq_along))
  x
}
for (t in c("distributions", "log_lik")) {
  b <- rd(sprintf("dev/cdfs-out/bin2-rec/base-rec-frmtmb-%s.tsv", t))
  n <- rd(sprintf("dev/cdfs-out/bin2-rec/rec-frmtmb-%s.tsv", t))
  m <- merge(b[, c("key", "V5", "V7")], n[, c("key", "V5", "V7")], by = "key",
             all = TRUE, suffixes = c(".base", ".lane"))
  ch <- m[is.na(m$V5.base) | is.na(m$V5.lane) | m$V5.base != m$V5.lane, ]
  cat(sprintf("== %s: %d records, %d moved\n", t, nrow(m), nrow(ch)))
  for (i in seq_len(nrow(ch))) {
    cat(sprintf("%-28s held %s -> %s | lane: %s\n", ch$key[i], ch$V5.base[i],
                ch$V5.lane[i], substr(ch$V7.lane[i], 1, 90)))
  }
}
