## Shared Lotka-Volterra model used by the RTMBode sampler probes.
## frmtmb-free: bare RTMB + RTMBode, so the reproduction can be filed upstream.
##
## Mirrors what frmtmb builds for a two-group ODE fit: log-scale rate
## parameters, one RTMBode::ode() call per group, lognormal observations.

suppressPackageStartupMessages({
    library(RTMB)
    library(RTMBode)
})

lv_dyn <- function(t, y, p) {
    ing <- p[1] * y[1] * y[2]
    list(c(p[2] * y[1] * (1 - y[1] / p[5]) - ing,
           ing * p[4] - p[3] * y[2]))
}

lv_truth <- c(rIng = 0.2, rGrow = 1.0, rMort = 0.2, assEff = 0.5, K = 10)
lv_times <- seq(0, 30, by = 1)
lv_inits <- list(c(Prey = 1, Predator = 2), c(Prey = 3, Predator = 1))

## Simulate `ngroup` independent trajectories with lognormal noise.
lv_data <- function(ngroup = 1L, sd = 0.1, seed = 1L) {
    set.seed(seed)
    lapply(seq_len(ngroup), function(g) {
        sol <- deSolve::ode(lv_inits[[g]], lv_times, lv_dyn, lv_truth)[, -1]
        sol * exp(rnorm(length(sol), sd = sd))
    })
}

## Objective: one RTMBode::ode() adjoint node per group.
lv_objective <- function(obs, silent = TRUE) {
    ngroup <- length(obs)
    likfun <- function(pl) {
        getAll(pl)
        p <- exp(logpar)
        nll <- 0
        for (g in seq_len(ngroup)) {
            sol <- RTMBode::ode(lv_inits[[g]], lv_times, lv_dyn, p,
                                atol = 1e-8, rtol = 1e-8)[, -1]
            nll <- nll - sum(dnorm(log(obs[[g]]), log(sol), exp(logsd), log = TRUE))
        }
        nll
    }
    pars <- list(logpar = log(lv_truth), logsd = log(0.1))
    MakeADFun(likfun, pars, silent = silent)
}
