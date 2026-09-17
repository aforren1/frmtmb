.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_lp_basis
### Title: The design of a linear predictor over the coefficient vector
### Aliases: frm_lp_basis

### ** Examples

set.seed(1)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
fit <- frm(bf(y ~ x + (1 | g)), data = dd)
nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
str(lb$A)
lb$coef_names

# the covariance of the WHOLE grid, which predict() reduces to its
# diagonal
Sigma <- lb$A %*% lb$V %*% t(lb$A)
all.equal(sqrt(diag(Sigma)),
          predict(fit, newdata = nd, re_formula = NA, se.fit = TRUE)$se.fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
