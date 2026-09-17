.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.ode))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_lincmt
### Title: Analytic one- to three-compartment pharmacokinetics
### Aliases: frm_lincmt

### ** Examples

# One-compartment oral pharmacokinetics, 100 into the depot every 12
# hours, and the same schedule already at steady state.
doses <- data.frame(time = 0, state = "depot", value = 100,
                    ii = 12, addl = 3L)
frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
           times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
           events = doses)

frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
           times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
           events = data.frame(time = 0, state = "depot", value = 100,
                               ii = 12, ss = TRUE))

# The flip-flop case, where the textbook form returns NaN
frm_lincmt(parms = list(ka = 0.2, ke = 0.2, V = 10), times = 1:4,
           ncmt = 1, depot = TRUE, init = list(depot = 100))

# In a formula, with between-subject variability on both rates
set.seed(2026)
tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
dd <- data.frame(id = factor(rep(1:6, each = length(tt))),
                 time = rep(tt, 6), dose = 100)
ka <- exp(rnorm(6, 0, 0.3))[as.integer(dd$id)]
ke <- exp(rnorm(6, log(0.2), 0.25))[as.integer(dd$id)]
dd$conc <- 100 * ka / (10 * (ka - ke)) *
  (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)
## No test: 
fit <- frm(
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, init = list(depot = dose)),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
    gaussian(),
  data = dd, start = list(beta = c(0, log(0.25), log(8))))
fixef(fit)
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
