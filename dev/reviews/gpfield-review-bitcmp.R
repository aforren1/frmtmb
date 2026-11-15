# Review of lane gpfield: compare two bitwise RDS (gpfield-review-bitwise.R)
a <- commandArgs(TRUE)
o <- "C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/reviews/gpfield-review-out/"
x0 <- readRDS(paste0(o, "bit-", a[1], ".rds"))
y0 <- readRDS(paste0(o, "bit-", a[2], ".rds"))
for (n in names(x0)) {
  x <- x0[[n]]; y <- y0[[n]]
  cat(sprintf("BITCMP %-6s %s vs %s | est identical %s max|d| %.2e | ll identical %s | se identical %s max|dse| %.2e | warn identical %s\n",
              n, a[1], a[2], identical(x$est, y$est),
              max(abs(unlist(x$est) - unlist(y$est))), identical(x$ll, y$ll),
              identical(x$se, y$se), max(abs(x$se - y$se), na.rm = TRUE),
              identical(x$warn, y$warn)))
}
