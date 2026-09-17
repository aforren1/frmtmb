.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.learn))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: rlddm
### Title: Reinforcement learning with a drift-diffusion choice rule
### Aliases: rlddm

### ** Examples

d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 40,
                     seed = 1)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5), seed = 1)
fit <- frmtmb::frm(
  frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
             drift ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
  family = rlddm(subject = id, trial = trial), data = s[[1]])
frmtmb::fixef(fit)
head(frm_value_trace(fit))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
