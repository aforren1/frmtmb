# Bin 2 of the ported brms suite: the 0.77.0 ledger (HEAD,
# dev/rel078-log/bin2-ledger-077.tsv) against the 0.78.0
# consolidation's (dev/brmsport-bin2-ledger.tsv): outcome counts, every
# row whose outcome or class moved, every row only in one ledger, and
# every value check that moved. Adapted from dev/rel077-bin2-diff.R.
#   Rscript dev/rel078-bin2-diff.R > dev/rel078-log/bin2-diff.txt
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
b <- rd("dev/brmsport-bin2-ledger.tsv")
a <- rd("dev/rel078-log/bin2-ledger-077.tsv")
key <- function(x) paste(x$id, x$line)
stopifnot(!anyDuplicated(key(a)), !anyDuplicated(key(b)))
cat("rows: 0.77.0", nrow(a), "| 0.78.0", nrow(b), "\n")
print(table(before = paste(a$outcome, a$class, sep = " / ")))
print(table(after = paste(b$outcome, b$class, sep = " / ")))
only_a <- setdiff(key(a), key(b)); only_b <- setdiff(key(b), key(a))
cat("\nrows only in 0.77.0:", length(only_a), "| only in 0.78.0:",
    length(only_b), "\n")
for (k in only_a) cat("- 0.77.0 only:", k, "\n")
for (k in only_b) cat("- 0.78.0 only:", k, "\n")
sh <- intersect(key(a), key(b))
aa <- a[match(sh, key(a)), ]; bb <- b[match(sh, key(b)), ]
mv <- which(aa$outcome != bb$outcome | aa$class != bb$class)
cat("\nshared rows whose outcome or class moved:", length(mv), "\n")
writeLines(sprintf("- %s (%s%s, n %s): %s / %s -> %s / %s", aa$id[mv],
                   aa$target[mv],
                   ifelse(nzchar(aa$line[mv]), paste0(", line ", aa$line[mv]),
                          ""), aa$n[mv], aa$outcome[mv], aa$class[mv],
                   bb$outcome[mv], bb$class[mv]))
m <- which(aa$value_check != bb$value_check |
             aa$value_verdict != bb$value_verdict |
             aa$value_class != bb$value_class)
cat("\nvalue checks whose outcome or verdict moved:", length(m), "\n")
for (i in m) {
  cat(sprintf("- %s: %s / %s / %s -> %s / %s / %s\n", aa$id[i],
              aa$value_check[i], aa$value_verdict[i], aa$value_class[i],
              bb$value_check[i], bb$value_verdict[i], bb$value_class[i]))
}
