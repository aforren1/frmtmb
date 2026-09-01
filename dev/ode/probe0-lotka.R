# Probe 0: reproduce RTMBode's own ?ode Lotka-Volterra example verbatim,
# to establish that the adjoint ODE backend works on this machine before
# any frmtmb involvement.
library(RTMB)
library(RTMBode)

LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion <- rIng * Prey * Predator
    GrowthPrey <- rGrow * Prey * (1 - Prey / K)
    MortPredator <- rMort * Predator
    dPrey <- GrowthPrey - Ingestion
    dPredator <- Ingestion * assEff - MortPredator
    return(list(c(dPrey, dPredator)))
  })
}
pars <- c(rIng = 0.2, rGrow = 1.0, rMort = 0.2, assEff = 0.5, K = 10)
yini <- c(Prey = 1, Predator = 2)
times <- seq(0, 200, by = 1)

set.seed(1)
obs <- deSolve::ode(func = LVmod, y = yini, parms = pars, times = times)[, -1]
obs <- obs + rnorm(length(obs), sd = 1)

likfun <- function(p) {
  getAll(p)
  obs <- OBS(obs)
  sol <- ode(func = LVmod, y = yini, parms = pars, times = times,
             atol = 1e-8, rtol = 1e-8)
  obs %~% dnorm(mean = sol[, -1], sd = sdobs)
}

p <- list(pars = pars * 1.5, yini = yini * 1.5, sdobs = 1.5)
obj <- MakeADFun(likfun, p, silent = TRUE)
tm <- system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
print(tm)
cat("convergence:", opt$convergence, " objective:", opt$objective, "\n")
sdr <- sdreport(obj)
print(sdr)
cat("\n--- truth vs estimate ---\n")
print(rbind(truth = c(pars, yini, sdobs = 1),
            est = c(as.list(sdr, "Est")$pars, as.list(sdr, "Est")$yini,
                    as.list(sdr, "Est")$sdobs)))
cat("PROBE0 OK\n")
