.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: bandit2arm_delta
### Title: Rescorla-Wagner delta learning on a two-armed bandit
### Aliases: bandit2arm_delta

### ** Examples

d <- frm_task_design("bandit2arm", n_subject = 8, n_trial = 40,
                     seed = 1)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.4, tau = 3), seed = 1)[[1]]$choice
fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1,
                              tau ~ 1),
                   family = bandit2arm_delta(subject = id,
                                             trial = trial), data = d)
frmtmb::fixef(fit)
head(frm_value_trace(fit))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
