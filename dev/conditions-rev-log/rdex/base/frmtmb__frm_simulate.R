.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_simulate
### Title: Simulate responses from a formula and parameters
### Aliases: frm_simulate

### ** Examples

# power analysis: simulate from a design with chosen parameters
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                     newparams = list(b_Intercept = 1, b_x = 0.5,
                                      sigma = 0.7,
                                      sd_g__Intercept = 0.5),
                     nsim = 3, seed = 1)
head(sims)
# the same thing on the internal scale
par_template(bf(y ~ x + (1 | g)) + gaussian(), dd)   # the layout
sims2 <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                      newparams = list(beta = c(1, 0.5),
                                       betad = log(0.7),
                                       theta = log(0.5)),
                      nsim = 3, seed = 1)
# prior-predictive draws
pp <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                   prior = set_prior("normal(0, 1)", class = "b") +
                     set_prior("normal(0, 2)", class = "Intercept") +
                     set_prior("exponential(1)", class = "sd") +
                     set_prior("exponential(1)", class = "sigma"),
                   nsim = 4, seed = 1)
head(attr(pp, "pars"))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
