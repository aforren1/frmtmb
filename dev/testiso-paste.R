# Replace a marker line in a file with another file's lines.
#   Rscript dev/testiso-paste.R <target> <marker> <source>
a <- commandArgs(TRUE)
x <- readLines(a[1])
k <- which(x == a[2])
stopifnot(length(k) == 1L)
writeLines(c(x[seq_len(k - 1L)], readLines(a[3]), x[-seq_len(k)]), a[1])
