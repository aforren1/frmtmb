# From dev/nanse-fire-sum.R's output: on how many firings the lane names
# more, the same, or fewer parameters than base's sdreport() left with a
# non-finite standard error.
#   Rscript dev/nanse-fire-count.R <fire-block.txt>
x <- readLines(commandArgs(TRUE)[1])
x <- x[grepl("SE non-finite", x, fixed = TRUE)]
m <- regmatches(x, gregexpr("[0-9]+ of [0-9]+", x))
a <- vapply(m, function(v) as.integer(sub(" of.*", "", v[1])), 0L)
b <- vapply(m, function(v) as.integer(sub(" of.*", "", v[2])), 0L)
cat("firings", length(x), "; lane names more than base left non-finite:",
    sum(a > b), "; the same number:", sum(a == b), "; fewer:", sum(a < b),
    "\n")
