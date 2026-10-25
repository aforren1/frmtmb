# Bin 2 of the ported brms suite: the 0.74.0 ledger (HEAD,
# dev/rel075-log/bin2-ledger-074.tsv) against the 0.75.0
# consolidation's (dev/brmsport-bin2-ledger.tsv): outcome counts, every
# row whose outcome or class moved, and every value check that moved.
#   Rscript dev/rel075-bin2-diff.R > dev/rel075-log/bin2-diff.txt
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
b <- rd("dev/brmsport-bin2-ledger.tsv")
a <- rd("dev/rel075-log/bin2-ledger-074.tsv")
key <- function(x) paste(x$id, x$line)
stopifnot(setequal(key(a), key(b)), !anyDuplicated(key(a)))
print(table(before = a$outcome))
print(table(after = b$outcome))
bb <- b[match(key(a), key(b)), ]
mv <- which(a$outcome != bb$outcome | a$class != bb$class)
cat("\nrows whose outcome or class moved:", length(mv), "\n")
writeLines(sprintf("- %s (%s%s, n %s): %s / %s -> %s / %s", a$id[mv],
                   a$target[mv],
                   ifelse(nzchar(a$line[mv]), paste0(", line ", a$line[mv]),
                          ""), a$n[mv], a$outcome[mv], a$class[mv],
                   bb$outcome[mv], bb$class[mv]))
m <- which(a$value_check != bb$value_check |
             a$value_verdict != bb$value_verdict |
             a$value_class != bb$value_class)
cat("\nvalue checks whose outcome or verdict moved:", length(m), "\n")
for (i in m) {
  cat(sprintf("- %s: %s / %s / %s -> %s / %s / %s\n", a$id[i],
              a$value_check[i], a$value_verdict[i], a$value_class[i],
              bb$value_check[i], bb$value_verdict[i], bb$value_class[i]))
}
