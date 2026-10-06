# What a gp(x, by = f) fit prints: summary(), VarCorr(), confint_varcorr()
# and the brms-shaped coef tables, as a reader sees them.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                f = factor(rep(c("a", "b", "c"), length.out = n)))
d$y <- 0.5 + ifelse(d$f == "a", sin(d$x), cos(d$x)) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x, by = f)), data = d)
print(summary(fit))
print(VarCorr(fit))
print(confint_varcorr(fit))
print(fit)
