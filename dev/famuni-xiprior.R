# Does the prior the scalar-xi note suggests apply? frm_sample() with
# set_prior("normal(0, 2.5)", class = "xi") on a GEV with a scalar xi.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261017); n <- 300
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + 0.8 * (rexp(n)^(-0.15) - 1) / 0.15
fit <- frm(y ~ x, family = gen_extreme_value(), data = d)
dr <- frm_sample(fit, prior = set_prior("normal(0, 2.5)", class = "xi"),
                 chains = 2, iter = 600, seed = 1, refresh = 0)
print(prior_summary(dr))
print(summary(dr))
