.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: anova.frmtmb_multiple
### Title: Pooled model comparison across imputations (D1, D2, D3)
### Aliases: anova.frmtmb_multiple

### ** Examples

set.seed(4)
n <- 60
imps <- lapply(1:4, function(i) {
  x <- rnorm(n)
  data.frame(y = rnorm(n, 1 + 0.5 * x), x = x, z = rnorm(n))
})
m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imps)
m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)
anova(m1, m0)
anova(m1, m0, method = "D1")



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
