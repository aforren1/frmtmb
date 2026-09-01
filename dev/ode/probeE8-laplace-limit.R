# Probe E8: minimal, frmtmb-free reproduction of the second-order limit.
# n independent decays dy_i/dt = -exp(mu + u_i) y_i with u ~ N(0, 1)
# integrated out by Laplace. First-order AD is fine at any n; the Laplace
# inner Hessian goes NaN above a small state count. One n per process:
# a NaN ODE tape poisons later MakeADFun objects in the same session.
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1])
library(RTMB)
library(RTMBode)
times <- c(0, 1, 2, 3)
set.seed(1)
obs <- exp(-0.3 * rep(times[-1], each = n)) + rnorm(n * 3, 0, 0.1)

f <- function(p) {
  getAll(p)
  sol <- RTMBode::ode(y = rep(1, n), times = times,
                      func = function(t, y, q) list(-q * y),
                      parms = exp(mu + u), method = "lsoda")
  -sum(dnorm(u, 0, 1, log = TRUE)) +
    sum((obs - as.vector(t(sol[-1, -1, drop = FALSE])))^2)
}
o <- MakeADFun(f, list(mu = log(0.3), u = numeric(n)), random = "u",
               silent = TRUE)
v <- suppressWarnings(tryCatch(o$fn(log(0.3)), error = function(e) NaN))
g <- suppressWarnings(tryCatch(o$gr(log(0.3)), error = function(e) NaN))
cat(sprintf("n_states = %2d : Laplace fn = %-12s gr = %s\n", n,
            format(v, digits = 8), format(g, digits = 6)))
