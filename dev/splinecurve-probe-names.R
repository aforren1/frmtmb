# Lane splinecurve, round 2: which public accessor names the grouping
# factor of a bar term and its fitted levels, so that a refusal can name
# the term whose unseen level carries variance. Exploratory.
lib <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4)
d <- data.frame(g = factor(rep(1:10, each = 20)), h = factor(rep(1:4, 50)),
                t = runif(200))
d$y <- sin(3 * d$t) + rnorm(10, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
fit <- frm(bf(y ~ s(t, k = 6) + (1 + t | g) + (1 | h) +
                s(t, g, bs = "fs", k = 4)),
           family = gaussian(), data = d)
str(ngrps(fit))
r <- ranef(fit)
str(lapply(r, dimnames), max.level = 2)
str(lapply(VarCorr(fit), names), max.level = 2)
