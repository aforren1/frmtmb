# Lane optrap: the objective differences, lane minus base, of the fits
# dev/optrap-fitdiff.R calls the same optimum, in named test files.
#   Rscript dev/optrap-samediff.R <fitdiff tsv> <pattern>
a <- commandArgs(TRUE)
M <- read.delim(a[1], quote = "")
x <- M[grepl(a[2], M$tag) & M$class == "same optimum", ]
x$rel <- x$d / pmax(1, abs(x$obj.b))
x <- x[order(-abs(x$rel)), c("tag", "obj.b", "d", "rel")]
print(utils::head(x, 20), digits = 6, row.names = FALSE)
cat("fits:", nrow(x), "| lane lower:", sum(x$d < 0), "| lane higher:",
    sum(x$d > 0), "| equal:", sum(x$d == 0), "\n")
