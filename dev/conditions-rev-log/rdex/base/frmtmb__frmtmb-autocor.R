.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-autocor
### Title: Within-group residual correlation (R-side autocorrelation)
### Aliases: frmtmb-autocor

### ** Examples

set.seed(1)
d <- expand.grid(week = 1:5, subj = factor(1:30))
d$x <- rnorm(150)
e <- as.vector(vapply(1:30, function(i) {
  as.vector(stats::filter(rnorm(5), 0.6, "recursive"))
}, numeric(5)))
d$y <- 1 + 0.5 * d$x + e

fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
           data = d)
summary(fit)
autocor_matrix(fit)

# compound symmetry, and the unstructured correlation over the five
# weeks
frm(bf(y ~ x + cosy(week, subj)) + gaussian(), data = d)
frm(bf(y ~ x + unstr(week, subj)) + gaussian(), data = d)

# a random intercept alongside the correlated residual is allowed
frm(bf(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE)) + gaussian(),
    data = d)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
