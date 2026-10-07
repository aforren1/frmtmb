# The ported brms ledger before (0.68.1) and after this round: outcome
# counts and every row whose outcome, class or reason moved.
#   Rscript dev/rel069-ledger-diff.R > dev/rel069-log/ledger-diff.txt
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
a <- rd("dev/rel069-log/ledger-before.tsv")
b <- rd("dev/brmsport-ledger.tsv")
stopifnot(identical(sort(a$id), sort(b$id)))
lv <- c("pass", "cannot transfer", "divergence", "defect")
cat("| outcome | 0.68.1 | 0.69.0 |\n|---|---|---|\n")
for (o in lv) cat("|", o, "|", sum(a$outcome == o), "|", sum(b$outcome == o),
                  "|\n")
cat("| **total** | **", nrow(a), "** | **", nrow(b), "** |\n\n", sep = "")
b <- b[match(a$id, b$id), ]
mv <- which(a$outcome != b$outcome | a$class != b$class | a$reason != b$reason)
cat("Rows whose outcome, class or reason moved:", length(mv), "\n\n")
for (i in mv) {
  cat(sprintf("- `%s`: %s / %s -> %s / %s\n", a$id[i], a$outcome[i],
              ifelse(nzchar(a$class[i]), a$class[i], "-"), b$outcome[i],
              ifelse(nzchar(b$class[i]), b$class[i], "-")))
}
