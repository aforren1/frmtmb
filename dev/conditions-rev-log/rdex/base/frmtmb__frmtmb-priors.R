.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-priors
### Title: Prior objects, addressed by internal parameter name
### Aliases: frmtmb-priors prior_normal prior_t prior_logistic prior_gamma
###   prior_inv_gamma prior_beta prior_lkj

### ** Examples

# the objects themselves are cheap descriptions
prior_normal(0, 2)
prior_t(df = 3, location = 0, scale = 1)

set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)

# names are internal parameter names, or whole components. theta_1
# is a log-SD, so a normal there is a lognormal on the SD.
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           prior = list(beta = prior_normal(0, 5),
                        theta_1 = prior_t(3, 0, 1)))
prior_summary(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
