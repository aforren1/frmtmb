# Lane wt-defects: the before and after of the ported brms ledger,
# generated so that no count in dev/defects-findings.md is typed.
#   Rscript dev/defects-ledger-diff.R > dev/defects-log/ledger-diff.md
# Before: dev/defects-ledger-before.tsv, the ledger recorded at 0.65.0
# against rellib-r3 (dev/defects-log-record-r3.txt). After:
# dev/brmsport-ledger.tsv as regenerated from the lane library.
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
b <- rd("dev/defects-ledger-before.tsv")
a <- rd("dev/brmsport-ledger.tsv")
stopifnot(identical(sort(b$id), sort(a$id)), nrow(a) == 494L)
m <- merge(b[c("id", "outcome", "class")], a[c("id", "outcome", "class")],
           by = "id", suffixes = c(".b", ".a"))
lv <- c("pass", "defect", "divergence", "cannot transfer")
tb <- table(factor(m$outcome.b, lv))
ta <- table(factor(m$outcome.a, lv))
cat("| outcome | 0.65.0 (rellib-r3) | lane build |\n|---|---|---|\n")
for (o in lv) cat(sprintf("| %s | %d | %d |\n", o, tb[[o]], ta[[o]]))
cat(sprintf("| **total** | **%d** | **%d** |\n\n", sum(tb), sum(ta)))
cat(sprintf("Bin 1 passes: %d of 494 before, %d of 494 after.\n\n",
            tb[["pass"]], ta[["pass"]]))
mv <- m[m$outcome.b != m$outcome.a | m$class.b != m$class.a, ]
mv <- mv[order(mv$outcome.b, mv$outcome.a, mv$id), ]
cat(sprintf("Rows whose outcome or class moved: %d.\n\n", nrow(mv)))
cat("| row | before | after |\n|---|---|---|\n")
lab <- function(o, cl) if (nzchar(cl)) paste0(o, " (", cl, ")") else o
for (i in seq_len(nrow(mv))) {
  cat(sprintf("| `%s` | %s | %s |\n", mv$id[i],
              lab(mv$outcome.b[i], mv$class.b[i]),
              lab(mv$outcome.a[i], mv$class.a[i])))
}
d <- a[a$outcome == "defect", ]
cat(sprintf("\nDefect rows left: %d.\n\n", nrow(d)))
cat("| row | class |\n|---|---|\n")
for (i in seq_len(nrow(d))) cat(sprintf("| `%s` | %s |\n", d$id[i], d$class[i]))
