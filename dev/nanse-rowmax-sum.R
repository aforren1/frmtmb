# The measurement behind se_empty_row: over every firing of the
# standard-error warning in a dev/nanse-diag.sh directory, the largest
# Hessian entry of each parameter's row (optimizer units), split by
# whether the 1e-6 floor emptied it.
#   Rscript dev/nanse-rowmax-sum.R <dir>
d <- commandArgs(TRUE)[1]
fs <- list.files(d, pattern = "fire$", full.names = TRUE)
x <- unlist(lapply(fs, readLines))
rm <- grep("^  rowmax: ", x, value = TRUE)
rm <- rm[!grepl("NA$", rm)]
tok <- unlist(strsplit(sub("^  rowmax: ", "", rm), " "))
tok <- tok[grepl("=", tok, fixed = TRUE)]
val <- as.numeric(sub("^.*=([^*]*)[*]?$", "\\1", tok))
lost <- grepl("[*]$", tok)
# a row the logged finite-difference Hessian has as NaN (the c0 ridge)
# was judged on the exact Hessian, which is not logged
cat("rows not finite in the logged Hessian:", sum(!is.finite(val)), "\n")
lost <- lost[is.finite(val)]
val <- val[is.finite(val)]
lo <- val <= 1e-6
cat("firings with a Hessian logged:", length(rm), " rows:", length(val), "\n")
cat("rows at or below 1e-6 (emptied):", sum(lo), "; largest of them:",
    format(max(val[lo]), digits = 3), "\n")
cat("rows above 1e-6:", sum(!lo), "; smallest of them:",
    format(min(val[!lo]), digits = 3), "\n")
cat("  of those, lost by another rule:", sum(!lo & lost), "; kept:",
    sum(!lo & !lost), "; smallest kept:",
    format(min(val[!lo & !lost]), digits = 3), "\n")
