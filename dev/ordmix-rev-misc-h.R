.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib", "C:/Users/adf44/source/r/rellib-r5", "C:/Users/adf44/AppData/Local/R/win-library/4.6")); suppressPackageStartupMessages(library(frmtmb)); set.seed(20261095); n <- 300; d <- data.frame(x = rnorm(n), z = rnorm(n), ce = sample(0:1, n, TRUE)); lat <- ifelse(rbinom(n, 1, 0.4) == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n); d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5); d$g <- d$z + rnorm(n)
## (h) the multivariate fit is the sum of its two univariate fits: the
## responses share no parameter and no residual correlation
fm <- suppressWarnings(frm(mvbf(bf(y ~ x) + mixture(cumulative(), cumulative()),
                                bf(g ~ z) + gaussian()), data = d))
f_y <- suppressWarnings(frm(bf(y ~ x), family = mixture(cumulative(),
                                                        cumulative()),
                            data = d))
f_g <- frm(bf(g ~ z), family = gaussian(), data = d)
cat(sprintf("(h) mv logLik %.10f; univariate sum %.10f; diff %.3g\n",
            as.numeric(logLik(fm)),
            as.numeric(logLik(f_y)) + as.numeric(logLik(f_g)),
            as.numeric(logLik(fm)) - as.numeric(logLik(f_y)) -
              as.numeric(logLik(f_g))))
cat("(h) mv fixef rows:", paste(rownames(fixef(fm)), collapse = " "), "\n")
P1 <- fitted(fm, resp = "y")
P2 <- fitted(f_y)
cat(sprintf("(h) fitted(resp = 'y') vs univariate: max|diff| %.3g\n",
            max(abs(P1[, "Estimate", ] - P2[, "Estimate", ]))))
