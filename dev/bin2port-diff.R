# The bin-2 rows lane bin2port moved: the ledger of rellib-r9 before the
# lane (dev/bin2port-before-log/ledger.tsv, written by the driver on the
# unchanged tree) against the ledger after (dev/brmsport-bin2-ledger.tsv),
# printed as markdown for dev/bin2port-findings.md.
#
#   Rscript dev/bin2port-diff.R > dev/bin2port-out/diff.md
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
a <- rd("dev/bin2port-before-log/ledger.tsv")
b <- rd("dev/brmsport-bin2-ledger.tsv")
stopifnot(setequal(a$id, b$id))
b <- b[match(a$id, b$id), ]
lab <- function(d) {
  ifelse(d$outcome == "pass", paste0("pass (", d$kind, ")"),
         paste0(d$outcome, " (", d$class, ")"))
}
la <- lab(a)
lb <- lab(b)
mv <- which(la != lb)
cat("Rows moved:", length(mv), "statements,", sum(as.integer(b$n[mv])),
    "assertions\n\n")
cat("| row | n | before | after |\n|---|---|---|---|\n")
for (i in mv) {
  cat(sprintf("| `%s` | %s | %s | %s |\n", a$id[i], a$n[i], la[i], lb[i]))
}
cat("\nMoves by kind (assertions):\n\n| before | after | assertions |\n",
    "|---|---|---|\n", sep = "")
key <- paste(la[mv], "->", lb[mv])
tab <- tapply(as.integer(b$n[mv]), key, sum)
for (k in names(tab)) {
  p <- strsplit(k, " -> ", fixed = TRUE)[[1]]
  cat(sprintf("| %s | %s | %d |\n", p[1], p[2], tab[[k]]))
}
vc <- b[nzchar(b$value_check), ]
cat("\nValue checks on shape passes:", nrow(vc), "statements;",
    paste(sprintf("%s %d", names(table(vc$value_verdict)),
                  table(vc$value_verdict)), collapse = ", "), "\n")
