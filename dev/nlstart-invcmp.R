# The review's scale sweep (dev/review-nlstart-fix2-inv.R) on this build:
# logLik at s = 1, 1e-4, 1e-6 against each other, against the previous
# lane build (the review's after log) and 0.79.1 at 1e-6
# (fix2-small2-before.log is the review's base at 1e-6).
# Usage: Rscript dev/nlstart-invcmp.R <new log>
a <- commandArgs(TRUE)
rd <- function(f) {
  x <- read.table(text = grep("^ROW ", readLines(f), value = TRUE),
                  stringsAsFactors = FALSE)
  setNames(x[, 2:10], c("body", "fam", "seed", "l1", "l4", "l6", "w1", "w4",
                        "w6"))
}
A <- rd(a[1])
P <- rd("dev/review-nlstart-out/fix2-inv-after.log")
cat(readLines(a[1], 1), "\n")
best <- pmax(A$l1, A$l4, A$l6, P$l1, P$l4, P$l6, na.rm = TRUE)
cat(sprintf("rows %d; 1e-4 vs 1 differ > 1e-4: %d; 1e-6 vs 1 differ > 1e-4: %d\n",
            nrow(A), sum(abs(A$l4 - A$l1) > 1e-4), sum(abs(A$l6 - A$l1) > 1e-4)))
for (k in c("l1", "l4", "l6")) {
  cat(sprintf("%s: below the best of every arm by > 1e-4: now %d, previous build %d; now lower than previous by > 1e-4: %d; warned now %d\n",
              k, sum(best - A[[k]] > 1e-4), sum(best - P[[k]] > 1e-4),
              sum(A[[k]] - P[[k]] < -1e-4),
              sum(A[[sub("l", "w", k)]] > 0)))
}
