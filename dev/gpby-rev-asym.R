# Reviewer: extra_cov's asymmetry when two different doubles share a
# 15-significant-digit position key (0.1 + 0.2 and 0.3 + 7 - 7 style).
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5); d <- data.frame(x = sort(runif(40, 0, 4)))
d$y <- sin(d$x) + rnorm(40, 0, 0.2)
fit <- frm(bf(y ~ gp(x)), data = d)
a <- 6.1 + 0.2; b <- 6.3
nd <- data.frame(x = c(a, b, 7))
cat("a == b:", a == b, "| keys:", frmtmb:::pos_rowkey(as.matrix(nd$x)), "\n")
E <- frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)$extra_cov
cat(sprintf("isSymmetric %s | max |E - t(E)| %.3e | max |E| %.3e\n",
            isSymmetric(E), max(abs(E - t(E))), max(abs(E))))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d2 <- gpby_data()
f2 <- frm(bf(y ~ gp(x)), data = d2)
nd2 <- data.frame(x = c(7, 7 + 1e-12, 7 + 1e-9, 7 + 1e-6, 0.1 + 0.2, 0.3,
                        8, 8))
cat("0.3 observed:", any(d2$x == 0.3), "| keys:",
    frmtmb:::pos_rowkey(as.matrix(nd2$x)), "\n")
E2 <- frm_lp_basis(f2, newdata = nd2, extra_cov = TRUE)$extra_cov
D <- abs(E2 - t(E2))
w <- which(D == max(D), arr.ind = TRUE)[1, ]
cat(sprintf("near grid: max |E - t(E)| %.3e at rows %d,%d | E there %.6e vs %.6e\n",
            max(D), w[1], w[2], E2[w[1], w[2]], E2[w[2], w[1]]))
