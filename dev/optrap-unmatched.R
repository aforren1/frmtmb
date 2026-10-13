# Lane optrap: where the fits dev/optrap-fitdiff.R could not match sit.
#   Rscript dev/optrap-unmatched.R <fitdiff tsv>
M <- read.delim(commandArgs(TRUE)[1], quote = "")
x <- M[M$class %in% c("base only", "lane only"), ]
print(sort(table(x$tag), decreasing = TRUE))
y <- M[grepl("test-se-check", M$tag) & M$class == "lane lower (better)", ]
print(y[, c("tag", "obj.b", "code.b", "obj.n", "code.n", "msg.n")])
