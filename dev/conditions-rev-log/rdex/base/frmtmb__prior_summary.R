.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: prior_summary
### Title: Priors used in a fit
### Aliases: prior_summary prior_summary.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# what was actually applied, after set_prior() was matched to the
# coefficients of this design
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           prior = set_prior("normal(0, 1)", class = "b") +
                    set_prior("exponential(1)", class = "sd"))
prior_summary(fit)

# a plain maximum-likelihood fit reports that it had none
prior_summary(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
