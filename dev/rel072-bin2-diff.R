# Bin 2 of the ported brms suite: the 0.71.0 ledger (HEAD,
# dev/rel072-log/bin2-ledger-071.tsv), lane bin2port's ledger on
# 0.71.0 (dev/rel072-log/bin2-ledger-lane.tsv) and the 0.72.0
# consolidation's (dev/brmsport-bin2-ledger.tsv): outcome counts and
# every row whose outcome or class moved, against each.
#   Rscript dev/rel072-bin2-diff.R > dev/rel072-log/bin2-diff.txt
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
b <- rd("dev/brmsport-bin2-ledger.tsv")
key <- function(x) paste(x$id, x$line)  # bin2port filled the wiener rows' target
for (f in c("dev/rel072-log/bin2-ledger-071.tsv",
            "dev/rel072-log/bin2-ledger-lane.tsv")) {
  a <- rd(f)
  cat("==", f, "against dev/brmsport-bin2-ledger.tsv\n")
  stopifnot(setequal(key(a), key(b)), !anyDuplicated(key(a)))
  print(table(before = a$outcome))
  print(table(after = b$outcome))
  bb <- b[match(key(a), key(b)), ]
  mv <- which(a$outcome != bb$outcome | a$class != bb$class)
  cat("\nrows whose outcome or class moved:", length(mv), "\n")
  out <- sprintf("- %s (%s%s): %s / %s -> %s / %s", a$id[mv], a$target[mv],
                 ifelse(nzchar(a$line[mv]), paste0(", line ", a$line[mv]), ""),
                 a$outcome[mv], a$class[mv], bb$outcome[mv], bb$class[mv])
  if (length(mv) <= 80) writeLines(out) else writeLines(c(head(out, 80), "..."))
  cat("\n")
}
# the value checks beside the shape passes, against the lane's ledger
a <- rd("dev/rel072-log/bin2-ledger-lane.tsv")
bb <- b[match(key(a), key(b)), ]
m <- which(a$value_check != bb$value_check |
             a$value_verdict != bb$value_verdict |
             a$value_class != bb$value_class)
cat("value checks whose outcome or verdict moved against the lane:",
    length(m), "\n")
for (i in m) {
  cat(sprintf("- %s: %s / %s / %s -> %s / %s / %s\n", a$id[i],
              a$value_check[i], a$value_verdict[i], a$value_class[i],
              bb$value_check[i], bb$value_verdict[i], bb$value_class[i]))
}
