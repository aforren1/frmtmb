.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: ps
### Title: A penalized spline whose value a nonlinear body consumes
### Aliases: ps

### ** Examples

set.seed(1)
n_id <- 40
d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
sh <- stats::rnorm(n_id, 0, 0.05)[d$id]
d$y <- 2 + sin(2 * pi * (d$t + sh)) + stats::rnorm(nrow(d), 0, 0.1)
fit <- frm(bf(y ~ lev + ps(t + shift, k = 8),
              lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE),
           data = d)
fixef(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
