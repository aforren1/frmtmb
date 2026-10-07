# Lane setier: paste a generated block into dev/setier-findings.md in
# place of a PENDING-<tag> line.
#   Rscript dev/setier-fill.R <tag> <file with the block>
a <- commandArgs(TRUE)
f <- "dev/setier-findings.md"
x <- readLines(f)
i <- which(x == paste0("PENDING-", a[1]))
stopifnot(length(i) == 1L)
b <- paste0("    ", readLines(a[2]))
writeLines(c(x[seq_len(i - 1L)], b, x[-seq_len(i)]), f)
