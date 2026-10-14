# Reviewer of lane illse: whole-matrix readers outside vcov(): insight's
# get_varcov() on a lossy root fit (degree-8 raw polynomial), fresh fit.
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
f <- suppressWarnings(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = d))
cat("lossy", isTRUE(attr(ns$sdr_of(f)$cov_fixed_root, "lossy")), "\n")
r <- cap(insight::get_varcov(f))
cat("insight::get_varcov finite", all(is.finite(r$v)), "| warnings:", length(r$w), wtag(r$w), "\n")
lmf <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = d)
a <- c(1, 1.5^(1:8))
Vl <- vcov(lmf) * 191 / 200
cat("a' V a from get_varcov over lm's QR value:",
    sqrt(drop(t(a) %*% r$v[1:9, 1:9] %*% a)) / (predict(lmf, data.frame(x1 = 1.5, x2 = 1.5^2, x3 = 1.5^3, x4 = 1.5^4, x5 = 1.5^5, x6 = 1.5^6, x7 = 1.5^7, x8 = 1.5^8), se.fit = TRUE)$se.fit * sqrt(191 / 200)), "\n")
r2 <- cap(vcov_cluster(f, ~ x1 > 1.5))
