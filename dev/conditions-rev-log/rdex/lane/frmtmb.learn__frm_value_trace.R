.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_value_trace
### Title: Per-trial value estimates, prediction errors and choice
###   probabilities
### Aliases: frm_value_trace

### ** Examples

d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 30,
                     seed = 3)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.4, tau = 3), seed = 3)[[1]]$choice
fit <- frmtmb::frm(frmtmb::bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
                   family = bandit2arm_delta(subject = id,
                                             trial = trial), data = d)
tr <- frm_value_trace(fit)
head(tr)
# the identity that ties the trace to the likelihood
c(from_trace = sum(log(tr$p)), logLik = as.numeric(stats::logLik(fit)))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
