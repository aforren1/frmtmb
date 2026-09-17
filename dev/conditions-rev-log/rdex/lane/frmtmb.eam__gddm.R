.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: gddm
### Title: The generalized drift-diffusion family
### Aliases: gddm

### ** Examples

## No test: 
set.seed(1)
# a coarse grid, so that the example is quick; see gddm_control()
ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
dat <- gddm_simulate(400, mu = 2.5, bs = 3, ndt = 0.25, tau = 1,
                     bound = gddm_bound_exponential(), control = ctl)
dat$cond <- 1L
fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
           family = gddm(bound = gddm_bound_exponential(), control = ctl),
           data = dat)
fixef(fit)
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
