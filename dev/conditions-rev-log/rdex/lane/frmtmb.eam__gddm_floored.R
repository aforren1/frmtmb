.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: gddm_floored
### Title: How many rows the grid could not represent
### Aliases: gddm_floored

### ** Examples

## No test: 
set.seed(3)
ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
dat <- gddm_simulate(200, mu = 2, bs = 2.5, ndt = 0.25, control = ctl)
dat$cond <- 1L
fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
           family = gddm(control = ctl), data = dat)
gddm_floored(fit)
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
