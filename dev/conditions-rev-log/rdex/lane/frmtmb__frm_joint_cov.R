.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_joint_cov
### Title: The joint covariance of the fixed and random coefficients
### Aliases: frm_joint_cov

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
fit <- frm(bf(y ~ x + (1 | g)), data = dd)
jc <- frm_joint_cov(fit)
dim(jc$V)
table(jc$names)
head(jc$labels)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
