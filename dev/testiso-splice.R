# Replace one "## " section of dev/testiso-findings.md with a file's
# text, and paste a generated block into it.
#   Rscript dev/testiso-splice.R <heading prefix> <new section file>
#     [<marker> <file to paste at the marker>]
a <- commandArgs(TRUE)
f <- "dev/testiso-findings.md"
x <- readLines(f)
i <- which(startsWith(x, a[1]))
stopifnot(length(i) == 1L)
nxt <- which(startsWith(x, "## ") & seq_along(x) > i)
j <- if (length(nxt)) nxt[1] - 1L else length(x)
new <- readLines(a[2])
if (length(a) >= 4L) {
  k <- which(new == a[3])
  stopifnot(length(k) == 1L)
  new <- c(new[seq_len(k - 1L)], readLines(a[4]), new[-seq_len(k)])
}
writeLines(c(x[seq_len(i - 1L)], new, "", x[-seq_len(j)]), f)
