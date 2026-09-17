.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-links
### Title: Link functions
### Aliases: frmtmb-links

### ** Examples

set.seed(1)
d <- data.frame(x = rnorm(80))
d$y <- rbinom(80, 1, pnorm(0.4 + 0.8 * d$x))

# the mean on a probit rather than a logit
fixef(frm(bf(y ~ x), family = bernoulli(link = "probit"), data = d))$mu

# a link on a parameter that is not the mean
d$z <- rnorm(80, 1 + d$x, exp(0.2 + 0.3 * d$x))
frm(bf(z ~ x, sigma ~ x), family = student(link_sigma = "log"), data = d)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
