# Bin 2 of the ported brms suite before (lane remeasure's ledger on
# 0.70.0, dev/rel071-log/bin2-ledger-before.tsv) and after the 0.71.0
# consolidation (dev/brmsport-bin2-ledger.tsv): outcome counts and every
# row whose outcome or class moved.
#   Rscript dev/rel071-bin2-diff.R > dev/rel071-log/bin2-diff.txt
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
a <- rd("dev/rel071-log/bin2-ledger-before.tsv")
b <- rd("dev/brmsport-bin2-ledger.tsv")
key <- function(x) paste(x$id, x$target, x$line)
stopifnot(setequal(key(a), key(b)), !anyDuplicated(key(a)))
print(table(before = a$outcome))
print(table(after = b$outcome))
b <- b[match(key(a), key(b)), ]
mv <- which(a$outcome != b$outcome | a$class != b$class)
cat("\nrows whose outcome or class moved:", length(mv), "\n")
out <- sprintf("- %s (%s%s): %s / %s -> %s / %s", a$id[mv], a$target[mv],
               ifelse(nzchar(a$line[mv]), paste0(", line ", a$line[mv]), ""),
               a$outcome[mv], a$class[mv], b$outcome[mv], b$class[mv])
writeLines(out)
