.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: cluster_scores
### Title: Per-cluster score matrix
### Aliases: cluster_scores

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:25, 8)))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + rnorm(25, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)

S <- cluster_scores(fit, ~ g)
dim(S)
# the scores add up to the gradient at the optimum
max(abs(colSums(S)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
