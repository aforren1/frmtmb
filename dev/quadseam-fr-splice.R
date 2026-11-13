# One-off for the fix round: replace quad_settle() and
# quad_settle_rounds in R/fit.R by dev/quadseam-fr-settle.txt.
f <- "C:/Users/adf44/source/r/frmtmb-wt-quadseam/R/fit.R"
x <- readLines(f)
a <- which(x == "#' Re-anchor the Gauss-Kronrod rule at the estimate until it is")
b <- which(x == "quad_settle_rounds <- 10L")
stopifnot(length(a) == 1L, length(b) == 1L, b > a)
new <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/quadseam-fr-settle.txt")
x <- c(x[seq_len(a - 1L)], new, x[-seq_len(b)])
writeLines(x, f)
cat("spliced", length(new), "lines\n")
