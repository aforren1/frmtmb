.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: simulate.frmtmb_fit
### Title: Simulate responses from a frmtmb fit
### Aliases: simulate.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# one column per draw; the seed used is attached
sims <- simulate(fit, nsim = 5, seed = 42)
str(sims)
attr(sims, "seed")

# re_formula = NA redraws the group effects, which is the right choice
# for a parametric bootstrap over new groups
sims_m <- simulate(fit, nsim = 5, re_formula = NA, seed = 42)
apply(sims_m, 2, var) > apply(sims, 2, var)

# a posterior-predictive check by hand: does the fit reproduce the
# share of zeros in the data?
mean(dd$y == 0)
colMeans(simulate(fit, nsim = 20, seed = 1) == 0)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
