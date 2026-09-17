.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: mvbf
### Title: Combine formulas into a multivariate model
### Aliases: mvbf set_rescor

### ** Examples

set.seed(2)
n <- 160
dd <- data.frame(x = rnorm(n), g = factor(rep(1:16, 10)))
u <- cbind(rnorm(16, 0, 0.8), rnorm(16, 0, 0.8))
e <- rnorm(n)                      # a disturbance both responses see
dd$y1 <- 1 + 0.5 * dd$x + u[dd$g, 1] + e + rnorm(n, 0, 0.5)
dd$y2 <- 2 - 0.3 * dd$x + u[dd$g, 2] + e + rnorm(n, 0, 0.5)

# each response keeps its own formula and family
fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
fixef(fit)

# rescor estimates the correlation of the residuals
fit_rc <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
              data = dd)
rescor_matrix(fit_rc)

# set_rescor() turns it on after the fact, and `+` also combines bf()s
mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE)
bf(y1 ~ x) + bf(y2 ~ x)

# |ID| correlates the random effects of the two responses
fit_id <- frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
                gaussian(), data = dd)
VarCorr(fit_id)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
