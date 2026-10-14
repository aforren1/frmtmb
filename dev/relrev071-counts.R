# Release reviewer of 0.71.0: recount the headline numbers of
# dev/round-20261014.md from the ledgers themselves (bin 1, bin 2), for
# 0.71.0 (the tree) and 0.70.x (dev/rel071-log/ledger-before.tsv,
# bin2-ledger-before.tsv).
rd <- function(f) utils::read.delim(f, colClasses = "character", quote = "")
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/"
for (v in c("0.70.x", "0.71.0")) {
  b1 <- rd(paste0(root, if (v == "0.71.0") "brmsport-ledger.tsv" else
    "rel071-log/ledger-before.tsv"))
  b2 <- rd(paste0(root, if (v == "0.71.0") "brmsport-bin2-ledger.tsv" else
    "rel071-log/bin2-ledger-before.tsv"))
  cat("\n==", v, "\n")
  cat("bin 1:", sum(b1$outcome == "pass"), "of", nrow(b1), "\n")
  print(table(b1$outcome))
  b2$w <- as.integer(b2$n)
  for (bn in c("2T", "2P")) {
    x <- b2[b2$bin == bn, ]
    cat("bin", bn, ":", sum(x$outcome == "pass"), "of", nrow(x), "\n")
  }
  print(xtabs(w ~ outcome + ct_kind, b2))
  dec1 <- sum(b1$outcome == "cannot transfer" &
                !(b1$class %in% c("absent", "fixture", "upstream",
                                  "adapter", "reference", "not-ported")))
  dec2 <- sum(b2$w[b2$outcome == "cannot transfer" & b2$ct_kind == "decision"])
  pass <- sum(b1$outcome == "pass") + sum(b2$w[b2$outcome == "pass"])
  tot <- nrow(b1) + sum(b2$w)
  cat(sprintf("combined %d of %d (%.1f%%)\n", pass, tot, 100 * pass / tot))
  cat("bin-1 cannot-transfer classes:\n")
  print(table(b1$class[b1$outcome == "cannot transfer"]))
  cat(sprintf("decisions: bin 1 %d (by class), bin 2 %d (ct_kind)\n", dec1,
              dec2))
  cat(sprintf("rate %d of %d (%.1f%%)\n", pass, tot - dec1 - dec2,
              100 * pass / (tot - dec1 - dec2)))
}
