# Does fitted(newdata) of one mo() row depend on how many rows the new
# data hold? test-mo-coding.R:236 asserts it does not, bitwise, and
# failed on ubuntu-26.04 (OpenBLAS 0.3.32) by one ulp in the k = 4 row.
#
#   Rscript dev/ci-moprobe.R <lib1;lib2;...>
# Prints, for matprod "default" (BLAS) and "internal" (R's own loops),
# whether the 3-row and 4-row predictions agree bitwise, and a direct
# BLAS probe of the same 2-column product at 3 and 4 rows.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cat("OPENBLAS_CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
# the data of test-mo-coding.R's mo_code_data()
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
f1 <- frm(y ~ mo(k), data = d)
for (mp in c("default", "internal")) {
  options(matprod = mp)
  p4 <- fitted(f1, newdata = data.frame(k = 1:4))
  p3 <- fitted(f1, newdata = data.frame(k = 2:4))
  dif <- p3 - p4[2:4, ]
  cat(sprintf(paste("matprod=%-8s identical=%s max|diff|/|value|=%.3g",
                    "ulps(k=4 Estimate)=%g\n"),
              mp, identical(p3, p4[2:4, ]),
              max(abs(dif) / abs(p4[2:4, ])),
              (p3[3, 1] - p4[4, 1]) / (2^floor(log2(abs(p4[4, 1]))) *
                                        .Machine$double.eps)))
}
# the same shape of product with BLAS: an intercept column and one
# column of values, 3 rows against 4
options(matprod = "default")
set.seed(1)
bad <- 0L
for (r in 1:2000) {
  b <- rnorm(2)
  v <- runif(4, 0, 5)
  X4 <- cbind(1, v)
  y4 <- drop(X4 %*% b)
  y3 <- drop(X4[2:4, , drop = FALSE] %*% b)
  if (!identical(y3, y4[2:4])) bad <- bad + 1L
}
cat("BLAS probe: 3-row and 4-row products differ in", bad, "of 2000\n")
