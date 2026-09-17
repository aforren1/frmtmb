.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: anova.frmtmb_fit
### Title: Likelihood-ratio tests between nested frmtmb fits
### Aliases: anova.frmtmb_fit

### ** Examples

set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.5))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)

m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
m1 <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
anova(m0, m1)

# dropping a variance component puts the null on the boundary, so
# this p-value is conservative by up to a factor of two
m2 <- frm(bf(y ~ x) + gaussian(), data = dd)
anova(m2, m0)

# REML fits compare only when the fixed-effect designs agree, which
# is the case for a variance-component test
r0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
r1 <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd, REML = TRUE)
anova(r0, r1)
# differing designs are refused; refit = TRUE compares ML fits instead
rz <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, REML = TRUE)
try(anova(r0, rz))
anova(r0, rz, refit = TRUE)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
