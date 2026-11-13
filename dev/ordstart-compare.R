# Lane ordstart: the 40-fit table before (rellib-r20, run 0 of
# dev/ordstart-starts.R's base output) and after (the lane build's
# default fit), with the kind of each answer: step (a component with
# sharpen > -0.1), runoff, collapsed or proper; the fit-end degenerate
# warning; the new step-start warning; and the logLik change.
#   Rscript dev/ordstart-compare.R before.tsv after.tsv [before after ...]
args <- commandArgs(trailingOnly = TRUE)
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE)
B <- do.call(rbind, lapply(args[c(TRUE, FALSE)], rd))
A <- do.call(rbind, lapply(args[c(FALSE, TRUE)], rd))
B <- B[B$run == 0, ]
A <- A[A$run == 0, ]
kind <- function(X) {
  st <- pmax(X$sharpen1, X$sharpen2, na.rm = TRUE) > -0.1
  ro <- pmax(X$near1, X$near2, na.rm = TRUE) > 50
  co <- X$collapsed1 | X$collapsed2
  ifelse(st, "step", ifelse(ro & co, "runoff+coll", ifelse(ro, "runoff",
         ifelse(co, "collapsed", "proper"))))
}
B$kind <- kind(B)
A$kind <- kind(A)
M <- merge(B[, c("family", "seed", "ll", "kind", "degen_warn", "elapsed")],
           A[, c("family", "seed", "ll", "kind", "degen_warn", "step_warn",
                 "starts_gain", "elapsed")],
           by = c("family", "seed"), suffixes = c("_before", "_after"))
M$dll <- M$ll_after - M$ll_before
op <- options(width = 200)
print(M[order(M$family, M$seed), ], digits = 7, row.names = FALSE)
mv <- abs(M$dll) > 1e-6
cat("\nfits", nrow(M), "; logLik moved by > 1e-6:", sum(mv),
    "; up:", sum(M$dll > 1e-6), "; down:", sum(M$dll < -1e-6), "\n")
cat("kind before -> after:\n")
print(table(before = M$kind_before, after = M$kind_after))
cat("degenerate warning before", sum(M$degen_warn_before), "after",
    sum(M$degen_warn_after), "; step-start warning after",
    sum(M$step_warn), "\n")
cat("proper before and down after:",
    sum(M$kind_before != "step" & M$dll < -1e-6), "\n")
cat("elapsed seconds, sum before", sum(M$elapsed_before), "after",
    sum(M$elapsed_after), "\n")
