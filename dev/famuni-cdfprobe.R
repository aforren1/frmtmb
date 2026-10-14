.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
g <- expand.grid(y = c(-30, -3, -0.5, 0, 0.7, 2, 8, 50), mu = c(-1, 0.5), sigma = c(0.3, 1, 4),
                 xi = c(-0.9, -0.4, -1e-7, 0, 1e-9, 1e-4, 0.01, 0.3, 1.2))
fam <- gen_extreme_value()
xi <- expm1(log1p(g$xi)); sg <- exp(log(g$sigma))
inside <- 1 + xi * (g$y - g$mu) / sg > 0
F <- fam$lcdf(g$y, list(mu = g$mu, sigma = sg, xi = xi), list())
print(cbind(g, F)[!inside & !is.finite(F), ][1:10, ])
print(sum(!inside & !is.finite(F)))
