.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: drop1.frmtmb_fit
### Title: Single-term deletions
### Aliases: drop1.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100), z = rnorm(100),
                 g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x * z + (1 | g)) + gaussian(), data = dd)

# marginality keeps the main effects out of scope while x:z is in it
drop1(fit)
drop1(fit, test = "Chisq")

# name the terms to override the default scope
drop1(fit, scope = ~ x + z, test = "Chisq")



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
