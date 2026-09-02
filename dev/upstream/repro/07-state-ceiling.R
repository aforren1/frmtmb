## The augmented-state ceiling.
##
## RTMBode differentiates by integrating an augmented system.  Its size is
##   order 0: nstate
##   order 1: nstate + nstate*(nstate+nparms)
##   order 2: nstate + nstate*(nstate+nparms) + nstate*(nstate+nparms)^2
## so a Laplace approximation (which needs order 2) integrates a system that
## grows as the CUBE of the model size.  The stiff Livermore integrators then
## allocate a dense Jacobian and a real work array of length
##   lrw = 22 + 9*neq + neq^2   (lsoda, full Jacobian)
## which is quadratic in that already-cubic number.
##
## One n per process: a failed ODE tape poisons later MakeADFun objects in the
## same session.
##
## Usage: Rscript 07-state-ceiling.R <n> <method> <lib>

args <- commandArgs(TRUE)
n <- as.integer(args[1])
method <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })

## nstate = n, nparms = n (one rate per state, all of them active).
aug <- function(order, nstate, nparms)
    sum(nstate * (nstate + nparms)^(0:order))
neq2 <- aug(2, n, n)
lrw <- 22 + 9 * neq2 + neq2^2

times <- c(0, 1, 2, 3)
set.seed(1)
obs <- exp(-0.3 * rep(times[-1], each = n)) + rnorm(n * 3, 0, 0.1)

f <- function(p) {
    getAll(p)
    sol <- RTMBode::ode(y = rep(1, n), times = times,
                        func = function(t, y, q) list(-q * y),
                        parms = exp(mu + u), method = method)
    -sum(dnorm(u, 0, 1, log = TRUE)) +
        sum((obs - as.vector(t(sol[-1, -1, drop = FALSE])))^2)
}

## First order only: no random effects, so the tape stops at order 1.
o1 <- MakeADFun(f, list(mu = log(0.3), u = numeric(n)), silent = TRUE)
g1 <- suppressWarnings(tryCatch(o1$gr(o1$par)[1], error = function(e) NaN))

## Laplace: needs order 2.
t0 <- proc.time()[3]
o2 <- MakeADFun(f, list(mu = log(0.3), u = numeric(n)), random = "u", silent = TRUE)
v <- suppressWarnings(tryCatch(o2$fn(log(0.3)), error = function(e) NaN))
g <- suppressWarnings(tryCatch(o2$gr(log(0.3))[1], error = function(e) NaN))
el <- proc.time()[3] - t0

cat(sprintf("%-6s n=%-3d neq(order2)=%-8d lrw=%-12.4g (%.3g GB)  order1_gr=%-12s Laplace fn=%-12s gr=%-12s %.1fs\n",
            method, n, neq2, lrw, lrw * 8 / 1024^3,
            format(g1, digits = 6), format(v, digits = 8),
            format(g, digits = 6), el))
