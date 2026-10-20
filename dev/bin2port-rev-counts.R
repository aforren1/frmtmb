# Reviewer's recount of lane bin2port's bin-2 ledger, before and after,
# straight from the tsv files (no code shared with the ledger script).
#   Rscript dev/bin2port-rev-counts.R <after ledger> <before ledger>
a <- commandArgs(trailingOnly = TRUE)
rd <- function(f) utils::read.delim(f, quote = "", stringsAsFactors = FALSE,
                                    colClasses = "character",
                                    na.strings = NULL)
for (f in a) {
  L <- rd(f)
  L$n <- as.integer(L$n)
  cat("==", f, "rows", nrow(L), "assertions", sum(L$n), "\n")
  print(tapply(L$n, list(L$outcome, L$bin), sum))
  cat("pass", sum(L$n[L$outcome == "pass"]),
      " pass value", sum(L$n[L$outcome == "pass" & L$kind == "value"]),
      " pass shape", sum(L$n[L$outcome == "pass" & L$kind == "shape"]), "\n")
  ct <- L$outcome == "cannot transfer"
  print(tapply(L$n[ct], L$ct_kind[ct], sum))
  print(tapply(L$n[ct], L$class[ct], sum))
  if ("value_check" %in% names(L)) {
    sp <- L$outcome == "pass" & L$kind == "shape"
    cat("shape passes", sum(L$n[sp]), " with value check",
        sum(L$n[sp & nzchar(L$value_check)]), "\n")
    print(table(L$value_verdict[sp & nzchar(L$value_check)]))
  }
  b1 <- rd("dev/brmsport-ledger.tsv")
  b1$n <- if ("n" %in% names(b1)) as.integer(b1$n) else 1L
  dec1 <- sum(b1$n[b1$outcome == "cannot transfer" &
                     b1$ct_kind == "decision"])
  p1 <- sum(b1$n[b1$outcome == "pass"])
  dec2 <- sum(L$n[ct & L$ct_kind == "decision"])
  p2 <- sum(L$n[L$outcome == "pass"])
  tot <- sum(b1$n) + sum(L$n)
  cat("bin1 total", sum(b1$n), "pass", p1, "decisions", dec1, "\n")
  cat(sprintf("bins 1+2: %d of %d (%.1f%%); excl decisions %d of %d (%.1f%%)\n",
              p1 + p2, tot, 100 * (p1 + p2) / tot, p1 + p2,
              tot - dec1 - dec2, 100 * (p1 + p2) / (tot - dec1 - dec2)))
  cat(sprintf("over 709: %d of 709 (%.1f%%)\n", p1 + p2,
              100 * (p1 + p2) / 709))
}
