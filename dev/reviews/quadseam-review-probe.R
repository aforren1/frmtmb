.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = factor(rep(1:10, each = 10)), x = rnorm(100))
d$y <- rgamma(100, shape = 5, rate = 5 / exp(1 + 0.3 * d$x + rnorm(10, 0, 0.5)[d$g]))
f <- frm(y ~ x + (1 | g), family = Gamma(link = "log"), data = d,
         quadrature = TRUE)
str(f$estimates)
print(f$opt$par)
print(logLik(f))
# base 3477 case: the summary of the sd
set.seed(3)
d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
d$y <- rbinom(120, 1, plogis(0.2 + 0.5 * d$x))
f <- frm(y ~ x + (1 | g), family = bernoulli(), data = d, quadrature = TRUE)
print(VarCorr(f)[["g"]]$sd)
print(summary(f))
