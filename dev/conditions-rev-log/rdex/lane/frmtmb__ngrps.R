.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ngrps
### Title: Number of levels per random-effect grouping factor
### Aliases: ngrps ngrps.frmtmb_fit

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(100),
                 g = factor(rep(1:10, 10)),
                 h = factor(rep(1:4, each = 25)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)

# one count per distinct grouping factor
ngrps(fit)
# the count that decides whether a variance component is trustworthy,
# and the unit influence() deletes when given `groups`
ngrps(fit)[["h"]]



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
