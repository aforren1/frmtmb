.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: mixture_mvn
### Title: Multivariate gaussian mixture family
### Aliases: mixture_mvn

### ** Examples

set.seed(1)
Y <- rbind(matrix(rnorm(60, 0), ncol = 2),
           matrix(rnorm(60, 4), ncol = 2))
dd <- data.frame(row = seq_len(nrow(Y)))
dd$Y <- Y
fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
fixef(fit)
head(mixture_probs(fit))
# a shared spherical covariance (mclust's EII, k-means-like)
frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = "EII"), data = dd)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
