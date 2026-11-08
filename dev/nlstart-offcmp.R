# Do the cost run's probe-off fits equal 0.79.1's default fits, row by row?
# Usage: Rscript dev/nlstart-offcmp.R <costdir> <scandir of rellib-r18>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(d, tag) unlist(lapply(list.files(d, "[.]txt$", full.names = TRUE),
  function(f) grep(paste0("^", tag, " "), readLines(f), value = TRUE)))
x <- read.table(text = sub("^COST ", "", rd(a[1], "COST")))
s <- read.table(text = sub("^ROW ", "", rd(a[2], "ROW")))
k1 <- paste(x$V1, x$V2, x$V3, x$V4); k2 <- paste(s$V1, s$V2, s$V3, s$V4)
d <- x$V5 - s$V5[match(k1, k2)]
cat("rows", length(d), " probe-off vs 0.79.1 default: |diff| > 1e-8:",
    sum(abs(d) > 1e-8), " max |diff|", format(max(abs(d)), digits = 3), "\n")
print(data.frame(key = k1, off = x$V5, base = s$V5[match(k1, k2)])[abs(d) > 1e-8, ],
      digits = 12, row.names = FALSE)
