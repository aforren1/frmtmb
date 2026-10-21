# Reviewer: brms-32 (truncated count mean) and brms-33 (inv_gaussian
# survivor sign), from brms 2.23.0's own code.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("frmtmb")
m <- brms:::posterior_epred_trunc_discrete("pois", list(lambda = matrix(3, 2, 2)),
                                            lb = matrix(2, 2, 2), ub = matrix(6, 2, 2))
k <- 2:6
cat(sprintf("brms posterior_epred_trunc_discrete(pois, 3, [2, 6]): %.10f\n", m[1]))
cat(sprintf("mean over [2, 6] (brms's likelihood window): %.10f\n",
            sum(k * dpois(k, 3)) / sum(dpois(k, 3))))
k <- 3:6
cat(sprintf("mean over (2, 6] = [3, 6]: %.10f\n", sum(k * dpois(k, 3)) / sum(dpois(k, 3))))
print(body(brms:::posterior_epred_trunc_discrete))
cat("== brms-33: inv_gaussian_lccdf as brms writes it, y = 2, mu = 2, shape = 4\n")
y <- 2; mu <- 2; sh <- 4
a <- sqrt(sh) / sqrt(y) * (y / mu - 1)
b <- sqrt(sh) / sqrt(y) * (y / mu + 1)
brms_lccdf <- log1p(-(pnorm(a) - exp(2 * sh / mu) * pnorm(-b)))
brms_lcdf <- log(pnorm(a) + exp(2 * sh / mu) * pnorm(-b))
cat(sprintf("brms lccdf %.6f | log1m(exp(brms lcdf)) %.6f | statmod %.6f | frmtmb %.6f\n",
            brms_lccdf, log1p(-exp(brms_lcdf)),
            statmod::pinvgauss(y, mu, sh, lower.tail = FALSE, log.p = TRUE),
            ns$inverse_gaussian_cdf_edges(y, mu, sh)$lS))
