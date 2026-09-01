# Probe E7: the stacked solve under Laplace only. One subject count per
# process (a failed ODE tape poisons every later MakeADFun in the same
# session - see the E5/E6 discrepancy), so this script takes n_id on the
# command line and is driven in a loop from the shell.
args <- commandArgs(trailingOnly = TRUE)
n_id <- as.integer(args[1])
library(RTMB)
library(RTMBode)
source("dev/ode/pk-common.R")
d <- sim_pk(n_id = n_id, seed = 100 + n_id)
gi <- as.integer(d$id); ng <- nlevels(d$id)
conc <- d$conc; time <- d$time; id <- d$id; dose <- d$dose

pk_dyn_stack <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  m <- length(y) / 2
  h <- 1:m
  A <- y[h]; C <- y[m + h]
  list(c(-p[h] * A, p[h] * A / p[2 * m + h] - p[m + h] * C))
}
pk_ode_stack <- function(ka, ke, V, time, id, dose) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  ids <- as.integer(factor(id)); ns <- max(ids)
  first <- match(seq_len(ns), ids)
  grid <- sort(unique(c(0, time)))
  sol <- RTMBode::ode(y = c(dose[first], rep(0, ns)), times = grid,
                      func = pk_dyn_stack,
                      parms = c(ka[first], ke[first], V[first]),
                      method = "lsoda", atol = 1e-8, rtol = 1e-8)
  out <- numeric(length(time))
  for (s in seq_len(ns)) {
    idx <- which(ids == s)
    out[idx] <- sol[match(time[idx], grid), 1L + ns + s]
  }
  out
}
f <- function(p) {
  getAll(p)
  mu <- pk_ode_stack(exp(lka + u_ka[gi]), exp(lke + u_ke[gi]),
                     exp(rep(lV, length(conc))), time, id, dose)
  -sum(dnorm(u_ka, 0, exp(lsd_ka), log = TRUE)) -
    sum(dnorm(u_ke, 0, exp(lsd_ke), log = TRUE)) -
    sum(dnorm(conc, mu, exp(lsigma), log = TRUE))
}
o <- MakeADFun(f, list(lka = 0, lke = log(0.25), lV = log(8), lsigma = 0,
                       lsd_ka = log(0.3), lsd_ke = log(0.3),
                       u_ka = numeric(ng), u_ke = numeric(ng)),
               random = c("u_ka", "u_ke"), silent = TRUE)
v <- suppressWarnings(tryCatch(o$fn(o$par), error = function(e) NaN))
g <- suppressWarnings(tryCatch(max(abs(o$gr(o$par))), error = function(e) NaN))
cat(sprintf("stacked Laplace, n_id=%2d (%d states): fn = %-12s max|gr| = %s\n",
            n_id, 2 * n_id, format(v, digits = 8), format(g, digits = 6)))
