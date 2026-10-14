# How rellib-r8 tables default priors (fit and sample routes) for a beta
# family's phi and a categorical response. Output: same name .log
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(2)
n <- 100
dd <- data.frame(x = rnorm(n), g = gl(10, 10))
dd$p <- plogis(rnorm(n))
dd$f <- factor(sample(c("a", "b", "c"), n, TRUE))
print(default_prior(bf(p ~ x + (1 | g)), data = dd, family = Beta()))
print(default_prior(bf(p ~ x + (1 | g)), data = dd, family = Beta(),
                    route = "sample"))
print(default_prior(bf(f ~ x + (1 | ID | g)), data = dd,
                    family = categorical(), route = "sample"))
dd$y1 <- rnorm(n); dd$y2 <- rnorm(n)
print(default_prior(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE),
                    data = dd, family = gaussian(), route = "sample"))
