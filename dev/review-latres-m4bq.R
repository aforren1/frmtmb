# The seeds of review-latres-m4b.R whose logLik sits far off the exact maximum.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4b2"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
x$over <- x$ll - x$x_ll
b <- x[is.finite(x$over) & x$over > 1, c("design", "seed", "ll", "x_ll", "over", "code", "nwarn", "ldisp", "x_ldisp", "ls", "x_ls", "ar", "x_ar", "span")]
print(b[order(b$design, b$seed), ], digits = 4, row.names = FALSE)
cat("\nper design: share with logLik > exact max + 1, of fits\n")
print(aggregate(over ~ design, x[is.finite(x$over), ], function(v) c(n = length(v), over1 = sum(v > 1), under1 = sum(v < -1), med = median(v))), digits = 3)
# Laplace reads true where sderr and ar are within what the exact MLE's own spread allows
