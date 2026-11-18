# Lane ordstarts: splice a generated log into dev/ordstarts-findings.md
# at a placeholder line, so numbers are pasted from logs, not typed.
#   Rscript dev/ordstarts-splice.R <placeholder> <log file> [append.md]
# With append.md, its content is appended to the findings first.
a <- commandArgs(trailingOnly = TRUE)
fn <- "dev/ordstarts-findings.md"
x <- readLines(fn)
if (length(a) > 2) x <- c(x, readLines(a[3]))
i <- which(x == a[1])
stopifnot(length(i) == 1L)
b <- readLines(a[2])
x <- c(x[seq_len(i - 1L)], b, x[-seq_len(i)])
writeLines(x, fn)
cat("spliced", length(b), "lines at", a[1], "\n")
