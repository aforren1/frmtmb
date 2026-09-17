.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: profile.frmtmb_fit
### Title: Likelihood profiles
### Aliases: profile.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# parameter names are the confint() row names, with or without the
# parentheses, or a one-to-one natural-scale alias
rownames(confint(fit))
pr <- profile(fit, "theta_1")
identical(profile(fit, "(Intercept)"), profile(fit, "Intercept"))
plot(pr)
# TMB's confint() reads the interval off the profile
confint(pr)

# several parameters at once return a named list
prs <- profile(fit, c("x", "theta_1"))
names(prs)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
