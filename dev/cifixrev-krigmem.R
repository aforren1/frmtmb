# Reviewer: peak memory of gp_krig_cov() on a 2000-row grid, for
# gp(x, by = <numeric>) and for a plain gp(x) grid with every position
# repeated once. Usage: Rscript dev/cifixrev-krigmem.R <lib|base>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(17)
xg <- round(runif(120, 0, 10), 1)
dg <- data.frame(y = sin(xg) + rnorm(120, 0, 0.3), x = xg,
                 w = runif(120, 0.5, 2))
fw <- suppressWarnings(frm(bf(y ~ gp(x, by = w)), data = dg))
fg <- suppressWarnings(frm(bf(y ~ gp(x)), data = dg))
peak <- function(expr) {
  invisible(gc(reset = TRUE))
  force(expr)
  g <- gc()
  g[2, 6]  # max used, Vcells, in Mb
}
nd_w <- data.frame(x = seq(0.01, 9.99, length.out = 2000),
                   w = seq(0.6, 1.9, length.out = 2000))
nd_d <- data.frame(x = rep(seq(0.01, 9.99, length.out = 1000), 2))
nd_u <- data.frame(x = seq(0.01, 9.99, length.out = 2000))
for (k in 1:2) {
  cat("by = numeric, 2000 distinct: peak Vcells Mb",
      peak(frm_lp_basis(fw, newdata = nd_w, extra_cov = TRUE)), "\n")
  cat("plain, 2000 rows at 1000 positions: peak Vcells Mb",
      peak(frm_lp_basis(fg, newdata = nd_d, extra_cov = TRUE)), "\n")
  cat("plain, 2000 distinct: peak Vcells Mb",
      peak(frm_lp_basis(fg, newdata = nd_u, extra_cov = TRUE)), "\n")
}
