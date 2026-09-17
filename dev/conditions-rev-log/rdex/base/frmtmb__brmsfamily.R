.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: brmsfamily
### Title: Build a family by name, with a link for any of its parameters
### Aliases: brmsfamily

### ** Examples

# sigma on a softplus: gaussian() from 'stats' has no link_sigma
brmsfamily("gaussian", link_sigma = "softplus")

# a mean link stats::poisson() refuses, and brms's short names
brmsfamily("poisson", softplus)$link
brmsfamily("zi_poisson")$link_zi

set.seed(2)
d <- data.frame(x = rnorm(60))
d$y <- rnorm(60, 1 + d$x, exp(0.3 * d$x))
frm(bf(y ~ x, sigma ~ x),
    family = brmsfamily("gaussian", link_sigma = "softplus"), data = d)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
