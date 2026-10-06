# Replace one test_that() block of a test file by the text in a file.
# Usage: Rscript dev/gpby-splice-test.R <test file> <first line of the
# block to replace> <first line of the block after it> <new text file>
a <- commandArgs(TRUE)
s <- readLines(a[1])
i <- which(s == a[2])
j <- which(s == a[3])
stopifnot(length(i) == 1, length(j) == 1, i < j)
out <- c(s[seq_len(i - 1L)], readLines(a[4]), s[j:length(s)])
writeLines(out, a[1])
cat("replaced lines", i, "to", j - 1L, "\n")
