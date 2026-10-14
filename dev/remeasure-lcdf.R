# R8: a frmtmb family's $lcdf slot returns the CDF (R/families.R:113),
# while row_lpdf() of a censored row gives the right log CDF and log
# survivor values. inverse.gaussian, mu 1, shape 3.
#
#   Rscript dev/remeasure-lcdf.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
fam <- brmsfamily("inverse.gaussian")
x <- c(1.2, 3)
dp <- list(mu = c(1, 1), shape = c(3, 3))
p <- brms::pinv_gaussian(x, 1, 3)
cat("lcdf:", fam$lcdf(x, dp, list()), " pinv:", p, " log(pinv):", log(p),
    "\n")
cat("lccdf:", if (is.function(fam$lccdf)) fam$lccdf(x, dp, list()) else
  "none", "\n")
for (cc in c(-1, 1)) {
  cat("row_lpdf cens", cc, ":",
      row_lpdf(fam, x, x, dp, list(cens = c(cc, cc)), NULL), " expected:",
      if (cc == 1) log(1 - p) else log(p), "\n")
}
