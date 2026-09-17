.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_control
### Title: Control parameters for frmtmb fits
### Aliases: frmtmb_control

### ** Examples

set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# another optimizer, with its own control list
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           control = frmtmb_control(optimizer = "optim",
                                    optCtrl = list(maxit = 500)))
fit$opt$convergence

# a tighter gradient criterion, with restarts from the current
# optimum until it is met
frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
    control = frmtmb_control(grad_tol = 1e-4, restarts = 3))

# badly scaled predictors: fit an internally standardized copy first,
# then warm-start the reported fit from it
dd$xbig <- dd$x * 1e5
frm(bf(y ~ xbig + (1 | g)) + gaussian(), data = dd,
    control = frmtmb_control(autoscale = TRUE))

# the object is a plain list, so it can be built once and reused
ctrl <- frmtmb_control(check_nlev_1 = "ignore")
ctrl$optimizer



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
