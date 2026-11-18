# Lane drawdir: print chosen columns of probe tsvs for ids matching a regex.
# Usage: Rscript dev/drawdir-col.R <regex> <cols comma> <tsv> [<tsv> ...]
a <- commandArgs(TRUE)
d <- do.call(rbind, lapply(a[-(1:2)], function(f) {
  utils::read.delim(f, stringsAsFactors = FALSE)[, c("id", strsplit(a[2], ",")[[1]])]
}))
d <- d[grepl(a[1], d$id), ]
num <- vapply(d, is.numeric, TRUE)
d[num] <- lapply(d[num], signif, 4)
print(d, row.names = FALSE, right = FALSE)
