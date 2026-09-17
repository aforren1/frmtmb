.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.ode))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frm_ode
### Title: Solve an ODE once per group inside a nonlinear predictor
### Aliases: frm_ode

### ** Examples

# One-compartment oral pharmacokinetics with between-subject
# variability on the absorption and elimination rates.
#   dA/dt = -ka A            A(0) = dose
#   dC/dt =  ka A / V - ke C C(0) = 0
pk_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

set.seed(2026)
tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
n_id <- 6
dd <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                 time = rep(tt, n_id), dose = 100)
ka <- exp(rnorm(n_id, 0, 0.3))[as.integer(dd$id)]
ke <- exp(rnorm(n_id, log(0.2), 0.25))[as.integer(dd$id)]
dd$conc <- 100 * ka / (10 * (ka - ke)) *
  (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)

if (requireNamespace("RTMBode", quietly = TRUE)) {
  ## No test: 
  fit <- frm(
    bf(conc ~ frm_ode(pk_dyn,
                      init   = list(dose, 0),
                      times  = time,
                      parms  = list(exp(lka), exp(lke), exp(lV)),
                      group  = id,
                      states = c("depot", "central"),
                      output = "central"),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = dd, start = list(beta = c(0, log(0.25), log(8))))
  fixef(fit)
  
## End(No test)

  # Repeated dosing: 100 into the depot every 12 hours. The dose at
  # time 0 is the initial condition, the rest are events.
  doses <- data.frame(time = c(12, 24, 36), state = "depot",
                      value = 100)
  frm_ode(pk_dyn, init = list(100, 0), times = c(6, 18, 30, 42),
          parms = list(1, 0.2, 10), states = c("depot", "central"),
          output = "central", events = doses)

  # The same schedule written compactly, and already at steady state
  frm_ode(pk_dyn, init = list(0, 0), times = c(6, 18, 30, 42),
          parms = list(1, 0.2, 10), states = c("depot", "central"),
          output = "central",
          events = data.frame(time = 0, state = "depot", value = 100,
                              ii = 12, ss = TRUE))

  # An elimination rate that doubles after hour 6, carried forward
  # from the row it appears on. `tv` values follow `parms`, so this
  # dynamics reads ka at p[1], V at p[2] and the time-varying ke at
  # p[3].
  pk_tv <- function(t, y, p) {
    list(c(-p[1] * y[1], p[1] * y[1] / p[2] - p[3] * y[2]))
  }
  tt <- c(1, 3, 6, 9, 12)
  frm_ode(pk_tv, init = list(100, 0), times = tt, parms = list(1, 10),
          tv = list(ifelse(tt < 6, 0.2, 0.4)),
          states = c("depot", "central"), output = "central")
}



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
