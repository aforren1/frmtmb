.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb.learn-package
### Title: frmtmb.learn: Reinforcement-Learning Families for 'frmtmb'
###   Models
### Aliases: frmtmb.learn frmtmb.learn-package
### Keywords: internal

### ** Examples

# a reversal task, and the change in learning rate with an interval
d <- frm_task_design("reversal", n_subject = 10, n_trial = 40, seed = 5)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.35, tau = 3), seed = 5)[[1]]$choice
fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal, tau ~ 1),
  family = bandit2arm_delta(subject = id, trial = trial), data = d)
frmtmb::fixef(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
