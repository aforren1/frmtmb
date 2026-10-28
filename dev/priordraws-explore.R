# Layout of a frmtmb_draws object on rellib-r14: columns, stored
# pieces, the stanfit's warmup record.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] + rnorm(60)
fit <- frm(bf(y ~ x + z + (1 + x | g)), data = d, family = gaussian())
ds <- suppressMessages(frm_sample(fit, chains = 2, iter = 300, seed = 1,
                                  refresh = 0,
                                  prior = set_prior("normal(0, 1)")))
print(names(ds))
print(colnames(ds$draws))
print(dim(ds$draws))
print(ds$stanfit@sim[c("iter", "warmup", "warmup2", "n_save", "chains")])
print(str(ds$fit$prior))
print(prior_summary(ds))
pl <- ds$fit$prior
ri <- frmtmb::resolve_prior_input(ds$fit, pl)
str(ri$entries, max.level = 3)
