.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.sample))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_sample
### Title: Sample a model with NUTS
### Aliases: frm_sample

### ** Examples

## No test: 
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
summary(ds)
fixef(ds)
hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
           class = NULL)

# the same model sampled straight from the formula, with no ML fit.
# It reports the brms default priors it chose, and prior_summary()
# gives them back.
ds2 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
                  chains = 1, iter = 500, refresh = 0)
prior_summary(ds2)
fixef(ds2)

# a set_prior() specification takes over the classes it names and
# leaves the rest of the defaults alone
ds3 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
                  chains = 1, iter = 500, refresh = 0,
                  prior = set_prior("exponential(1)", class = "sd"))
prior_summary(ds3)
}
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
