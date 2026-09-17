.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: mixture_probs
### Title: Posterior class probabilities of a mixture fit
### Aliases: mixture_probs

### ** Examples

set.seed(3)
dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)

p <- mixture_probs(fit)
head(p)
rowSums(p)[1:3]                    # rows sum to one

# the hard assignment, and how well it recovers the truth
cl <- max.col(p)
table(cl, truth = rep(1:2, each = 80))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
