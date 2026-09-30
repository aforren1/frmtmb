# brms 2.23.0's own variables(), summary and default_prior for the
# equated mixture, seed 11, the data of dev/formula2-probe-equate-frm.R.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
fam <- mixture(gaussian, gaussian)
fit <- brm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d,
           chains = 1, iter = 400, refresh = 0, seed = 1)
cat("VARIABLES\n")
print(variables(fit))
print(summary(fit))
print(fixef(fit))
fit2 <- brm(bf(y ~ x, sigma1 = "sigma2", theta1 ~ 1), family = fam,
            data = d, chains = 1, iter = 400, refresh = 0, seed = 1)
cat("VARIABLES theta1 ~ 1\n")
print(variables(fit2))
