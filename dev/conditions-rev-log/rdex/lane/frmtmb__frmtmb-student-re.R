.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb-student-re
### Title: Student-t distributed random effects
### Aliases: frmtmb-student-re

### ** Examples

set.seed(1)
n <- 12
d <- data.frame(x = rnorm(20 * n), g = factor(rep(1:20, each = n)))
b <- rnorm(20)
b[20] <- b[20] + 6          # one outlying group
d$y <- 1 + 0.5 * d$x + b[d$g] + rnorm(20 * n)

fit_t <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))),
             family = gaussian(), data = d)
fit_n <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)

# the gaussian latent has to widen to cover the outlying group
VarCorr(fit_t)
VarCorr(fit_n)

# heavier tails, at the cost of a fixed nu
frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3))),
    family = gaussian(), data = d)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
