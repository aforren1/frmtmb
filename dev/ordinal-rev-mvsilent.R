# Reviewer: the loop-closure capture in draws_ordinal_cols(). Two
# equidistant ordinal responses with the same threshold count, cumulative
# (delta = exp(internal)) then sratio (delta = internal): is delta_y in the
# draws exp() of its internal column? Data seed 20261010, sampler seed 5.
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib", "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261010)
n <- 200
d <- data.frame(x = rnorm(n))
mk <- function() 1L + (stats::rlogis(n) + d$x > -1) + (stats::rlogis(n) + d$x > 0) +
  (stats::rlogis(n) + d$x > 1) + (stats::rlogis(n) + d$x > 2)
d$y <- mk(); d$y2 <- mk()
fit <- frm(bf(y ~ x) + bf(y2 ~ x), data = d,
           family = list(cumulative(threshold = "equidistant"),
                         sratio(threshold = "equidistant")))
cat("ML: summary delta rows\n"); print(summary(fit)$spec_pars)
ds <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1, iter = 200,
                                                   refresh = 0, seed = 5)))
M <- as.matrix(ds)
th <- M[, grep("^b_y_Intercept[[]", colnames(M))]
cat("draws: mean delta_y column", mean(M[, "delta_y"]),
    "; mean spacing of the stored b_y_Intercept columns", mean(th[, 2] - th[, 1]), "\n")
cat("draws: max |delta_y - (b_y_Intercept[2] - b_y_Intercept[1])|",
    max(abs(M[, "delta_y"] - (th[, 2] - th[, 1]))), "\n")
th2 <- M[, grep("^b_y2_Intercept[[]", colnames(M))]
cat("draws: max |delta_y2 - spacing|", max(abs(M[, "delta_y2"] - (th2[, 2] - th2[, 1]))), "\n")
