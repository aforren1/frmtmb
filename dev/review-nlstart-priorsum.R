# Review of lane nlstart: prior sweep, lane body against base body.
rd <- function(f) {
  cat("read", f, ":", readLines(f, 1), "\n")
  x <- read.table(text = grep("^ROW ", readLines(f), value = TRUE))
  setNames(x[, 2:10], c("case", "s", "psd", "seed", "ll", "conv", "nwarn",
                        "ll_lin", "conv_lin"))
}
A <- rd("dev/review-nlstart-out/prior-after.log")
B <- rd("dev/review-nlstart-out/prior-before.log")
d <- A$ll - B$ll
cat(sprintf("rows %d; lane body lower than base by > 1e-6: %d (range %.2g to %.2g), higher: %d\n",
            nrow(A), sum(d < -1e-6), min(-d[d < -1e-6]), max(-d[d < -1e-6]),
            sum(d > 1e-6)))
w <- d < -1e-6
print(data.frame(A[w, c("case", "s", "psd", "seed", "conv", "nwarn")],
                 base_conv = B$conv[w], diff = d[w]), row.names = FALSE)
