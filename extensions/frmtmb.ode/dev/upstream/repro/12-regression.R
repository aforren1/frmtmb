## Regression: the package's own documented example, and a numeric-path solve,
## must be unchanged by the patches.
##
## Usage: Rscript 12-regression.R <lib>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })
cat("RTMBode from:", dirname(system.file(package = "RTMBode")), "\n\n")

## ---- ?ode example, verbatim ----
LVmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        Ingestion <- rIng * Prey * Predator
        GrowthPrey <- rGrow * Prey * (1 - Prey/K)
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
tm <- system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))[3]
cat("== ?ode example ==\n")
cat("objective:", format(opt$objective, digits = 10), " convergence:", opt$convergence,
    " iterations:", opt$iterations, " seconds:", round(tm, 1), "\n")
sdr <- sdreport(obj)
est <- as.list(sdr, "Est")
cat("pars :", signif(est$pars, 6), "\n")
cat("yini :", signif(est$yini, 6), "\n")
cat("sdobs:", signif(est$sdobs, 6), "\n")
cat("max |gr| at optimum:", max(abs(obj$gr(opt$par))), "\n")

## ---- numeric path ----
cat("\n== numeric path (no advectors) ==\n")
s <- ode(y = yini, times = 0:10, func = LVmod, parms = pars)
cat("dim:", dim(s), " colnames:", colnames(s), "\n")
sref <- deSolve::ode(y = yini, times = 0:10, func = LVmod, parms = pars)
cat("max abs difference from deSolve::ode:", max(abs(s - sref)), "\n")

## ---- gradient against finite differences ----
cat("\n== gradient check ==\n")
g <- obj$gr(obj$par)
h <- 1e-6
fd <- sapply(seq_along(obj$par), function(i) {
    a <- b <- obj$par; a[i] <- a[i] + h; b[i] <- b[i] - h
    (obj$fn(a) - obj$fn(b)) / (2 * h)
})
cat("max relative error vs central differences:",
    signif(max(abs(as.numeric(g) - fd) / pmax(abs(fd), 1)), 4), "\n")
