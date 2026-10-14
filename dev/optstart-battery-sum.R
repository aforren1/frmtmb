# Lane optstart: the battery (dev/optstart-battery.R), before against
# after. A fit is short when its log likelihood is more than 1e-4 below
# the best of the reference package, the fit started at the true
# values and both builds' default fits. Per family and design: fits
# short (with the largest gap), at a nonzero code and warned, before and
# after, and the fits after that end lower than before.
#   Rscript dev/optstart-battery-sum.R <before label> <after label>
a <- commandArgs(TRUE)
rd <- function(lab) {
  fs <- Sys.glob(sprintf("dev/optstart-log/battery/%s-*.tsv", lab))
  do.call(rbind, lapply(fs, read.delim, stringsAsFactors = FALSE,
                        quote = ""))
}
B <- rd(a[1])
N <- rd(a[2])
M <- merge(B, N, by = c("family", "design", "seed"),
           suffixes = c(".b", ".n"))
cat("fits:", nrow(B), "before,", nrow(N), "after,", nrow(M), "matched\n")
# a log likelihood more than 50 above the fit started at the true
# values is a degenerate stop (huber() on 0.70.0 reached +1.4e5 on 400
# rows), not a maximum: counted on its own and kept out of the best
tr <- pmax(M$ll_truth.b, M$ll_truth.n, na.rm = TRUE)
M$deg.b <- !is.na(M$ll.b) & !is.na(tr) & M$ll.b > tr + 50
M$deg.n <- !is.na(M$ll.n) & !is.na(tr) & M$ll.n > tr + 50
best <- pmax(M$ll_ref.b, M$ll_truth.b, M$ll_truth.n,
             ifelse(M$deg.b, NA, M$ll.b), ifelse(M$deg.n, NA, M$ll.n),
             na.rm = TRUE)
M$gap.b <- M$ll.b - best
M$gap.n <- M$ll.n - best
M$refsrc <- ifelse(is.na(M$ll_ref.b), "truth start", "package")
s1 <- function(gap, code, nw) {
  sh <- !is.na(gap) & gap < -1e-4
  sprintf("%2d short%-10s %2d code %2d warned",
          sum(sh) + sum(is.na(gap)),
          if (any(sh)) sprintf(" (%.3g)", min(gap[sh])) else "",
          sum(is.na(code) | code != 0), sum(nw > 0))
}
cat(sprintf("\n%-13s %-9s %-4s | %-38s | %-38s\n", "family", "design", "n",
            "before (rellib-r8)", "after (lane)"))
for (fam in unique(M$family)) for (des in unique(M$design)) {
  x <- M[M$family == fam & M$design == des, ]
  if (!nrow(x)) next
  cat(sprintf("%-13s %-9s %-4d | %-38s | %-38s\n", fam, des, nrow(x),
              s1(x$gap.b, x$code.b, x$nwarn.b),
              s1(x$gap.n, x$code.n, x$nwarn.n)))
}
tot <- function(gap, code, nw) {
  c(short = sum(is.na(gap) | gap < -1e-4), code = sum(is.na(code) | code != 0),
    warned = sum(nw > 0))
}
cat("\ntotals before:", tot(M$gap.b, M$code.b, M$nwarn.b), "| after:",
    tot(M$gap.n, M$code.n, M$nwarn.n), "(short, nonzero code, warned)\n")
cat("degenerate stops (above the true-start fit by more than 50): before",
    sum(M$deg.b), "| after", sum(M$deg.n), "| errors: before",
    sum(is.na(M$ll.b)), "| after", sum(is.na(M$ll.n)), "\n")
cat("references: package", sum(M$refsrc == "package"), "| truth start only",
    sum(M$refsrc != "package"), "\n")
w <- M[!is.na(M$ll.n) & !is.na(M$ll.b) & M$ll.n < M$ll.b - 1e-4, ]
cat("\nafter lower than before by more than 1e-4:", nrow(w), "\n")
for (i in seq_len(nrow(w))) {
  cat(sprintf("  %s %s seed %d: before %.6f code %s | after %.6f code %s | best %.6f | %s\n",
              w$family[i], w$design[i], w$seed[i], w$ll.b[i], w$code.b[i],
              w$ll.n[i], w$code.n[i], best[match(TRUE,
              M$family == w$family[i] & M$design == w$design[i] &
                M$seed == w$seed[i])], substr(w$warn.n[i], 1, 60)))
}
u <- M[!is.na(M$gap.n) & M$gap.n < -1e-4, ]
cat("\nafter, short of the best:", nrow(u), "\n")
for (i in seq_len(nrow(u))) {
  cat(sprintf("  %s %s seed %d: gap %.4g code %s | %s\n", u$family[i],
              u$design[i], u$seed[i], u$gap.n[i], u$code.n[i],
              substr(u$warn.n[i], 1, 70)))
}
