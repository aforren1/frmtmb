.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: coef.frmtmb_fit
### Title: Per-group coefficients (fixed effects plus conditional modes)
### Aliases: coef.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# one row per group: the fixed effects with the modes added in
head(coef(fit)$g)
# which is fixef() plus ranef(), the lme4 identity
all.equal(coef(fit)$g[["(Intercept)"]],
          fixef(fit)$mu[["(Intercept)"]] + ranef(fit)$g[, 1],
          check.attributes = FALSE)

# without random effects there are no groups, so coef() is fixef()
coef(frm(bf(y ~ x) + gaussian(), data = dd))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
