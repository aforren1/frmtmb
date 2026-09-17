.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: update.frmtmb_fit
### Title: Update and refit a model
### Aliases: update.frmtmb_fit

### ** Examples

set.seed(3)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x), family = gaussian(), data = dd)

# a delta on the stored formula, in either spelling
fit2 <- update(fit, ~ . + z)
formula(fit2)
formula(update(fit, . ~ . + z))
fit3 <- update(fit, formula. = ~ . - x, newdata = dd[1:40, ])
nobs(fit3)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
