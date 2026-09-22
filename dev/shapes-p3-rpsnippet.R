.libPaths(c("C:/Users/adf44/source/r/shapes-lib","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
set.seed(3)
n <- 400
dd <- data.frame(trt = rbinom(n, 1, 0.5), centre = factor(rep(1:20, each = 20)))
u <- rnorm(20, 0, 0.3)
dd$t <- rweibull(n, shape = 1.4, scale = exp(1 - 0.5 * dd$trt + u[dd$centre]))
dd$censored <- as.integer(dd$t > 3); dd$t <- pmin(dd$t, 3)
fit <- frm(bf(t | cens(censored) ~ trt + (1 | c | centre),
              gamma1 ~ 1 + (1 | c | centre)),
           family = royston_parmar(df = 1), data = dd)
print(colnames(ranef(fit)[["centre"]]))
slope <- fixef_by_dpar(fit)$gamma1[["(Intercept)"]] +
  ranef(fit)[["centre"]][, "gamma1_Intercept"]
print(summary(slope))
print(rownames(ranef(fit)[["centre"]])[slope <= 0])
