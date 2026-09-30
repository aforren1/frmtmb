args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] == "before") character() else
  "C:/Users/adf44/source/r/wt-formula2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
print(fixef(fit))
print(variables(fit))
print(get_prior(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d))
str(fit$frame$linpreds[[3]][c("resp", "dpar", "idx", "par", "constant")])
print(names(fit$frame$linpreds))
print(fit$frame$par_template$betad)
print(logLik(fit))
