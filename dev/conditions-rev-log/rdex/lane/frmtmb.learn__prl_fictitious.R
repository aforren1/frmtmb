.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: prl_fictitious
### Title: Counterfactual (fictitious) updating on two options
### Aliases: prl_fictitious

### ** Examples

d <- frm_task_design("reversal", n_subject = 6, n_trial = 40, seed = 4)
# a reversal task pays -1 as often as +1
d$pay1 <- 2 * d$pay1 - 1
d$pay2 <- 2 * d$pay2 - 1
d$choice <- frm_task_simulate(
  prl_fictitious(subject = id, trial = trial), d,
  pars = list(alpha = 0.3, bias = 0, tau = 3), seed = 4)[[1]]$choice
fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ 1, bias ~ 1, tau ~ 1),
  family = prl_fictitious(subject = id, trial = trial), data = d)
frmtmb::fixef(fit)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
