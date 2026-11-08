# The review's small-covariate rows (dev/review-nlstart-fix-small.R) on this
# build, against the review's reference (best of 16 random starts and the
# truth, in its after log) and its rellib-r18 log.
# Usage: Rscript dev/nlstart-smallcmp.R <new log>
a <- commandArgs(TRUE)
rd <- function(f) {
  x <- read.table(text = grep("^ROW ", readLines(f), value = TRUE),
                  stringsAsFactors = FALSE, fill = TRUE)
  setNames(x[, 2:9], c("body", "fam", "seed", "ll", "conv", "nwarn", "multi",
                       "ref"))
}
A <- rd(a[1])
R <- rd("dev/review-nlstart-out/fix-small-after.log")
B <- rd("dev/review-nlstart-out/fix-small-before.log")
cat(readLines(a[1], 1), "\n")
best <- pmax(R$ref, R$ll, A$ll, B$ll, na.rm = TRUE)
A$short <- best - A$ll
B$short <- best - B$ll
for (fam in unique(A$fam)) for (body in unique(A$body)) {
  w <- A$fam == fam & A$body == body
  cat(sprintf("%-10s %-4s short > 1e-4: now %d (silent %d), base %d; lower than base: %d; shorts now: %s\n",
              fam, body, sum(A$short[w] > 1e-4), sum(A$short[w] > 1e-4 &
                A$nwarn[w] == 0), sum(B$short[w] > 1e-4),
              sum(A$ll[w] - B$ll[w] < -1e-4),
              paste(signif(A$short[w], 3), collapse = " ")))
}
