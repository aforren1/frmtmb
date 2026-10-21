# Reviewer, fix round 1: the zero_inflated_poisson window objective of
# cdfs-rev-win.R, row by row, along log mu
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
fam <- ns$as_frmtmb_family(zero_inflated_poisson())
y <- c(0, 0, 30, 2, 0)
av <- list(cens = c(2, -1, 2, 1, 0), cens_y2 = c(3, 0, 35, 2, 0))
f <- function(lm, lz = -1) ns$row_lpdf(fam, y, y, list(mu = exp(lm), zi = 1 / (1 + exp(-lz))), av, NULL)
for (lm in 1.5 + c(-1e-3, -1e-4, -1e-6, 0, 1e-6, 1e-4, 1e-3)) cat(sprintf("log mu %.7f: %s\n", lm, paste(format(f(lm), digits = 12), collapse = " ")))
cat("== row 3 (30, 35] alone over log mu, and its reference\n")
for (lm in c(1, 1.4, 1.5, 1.6, 2, 3)) {
  mu <- exp(lm)
  ref <- log(1 - 1 / (1 + exp(1))) + log(sum(dpois(31:35, mu)))
  cat(sprintf("log mu %.2f: row3 %.12g ref %.12g\n", lm, f(lm)[3], ref))
}
cat("== zi_poisson pair at q = 30 and 35, mu = exp(1.5)\n")
str(fam$cdf_log(c(30, 35), list(mu = exp(1.5), zi = 1 / (1 + exp(1))), list()))
