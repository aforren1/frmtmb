# Lane famend: the dpar names eval_dpars() gives a zero_inflated_poisson fit
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(31)
d <- data.frame(x = rnorm(100))
d$y <- ifelse(runif(100) < 0.2, 0L, rpois(100, exp(1 + 0.4 * d$x)))
f <- frm(y ~ x, family = zero_inflated_poisson(), data = d)
str(frmtmb:::eval_dpars(f)[[1]])
