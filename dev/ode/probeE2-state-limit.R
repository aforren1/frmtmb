# Probe E2: the stacked solve fails only under AD. Find the state-count
# threshold and whether it is integrator-specific. Minimal system, no
# frmtmb, no PK: n independent linear decays dy_i/dt = -r_i y_i, whose
# solution y_i(t) = exp(-r_i t) is known exactly.
library(RTMB)
library(RTMBode)

dyn <- function(t, y, p) { "c" <- RTMB::ADoverload("c"); list(-p * y) }
times <- seq(0, 5, length.out = 9)

run <- function(ns, method) {
  f <- function(par) {
    r <- exp(par$lr)
    sol <- RTMBode::ode(y = rep(1, ns), times = times, func = dyn, parms = r,
                        method = method, atol = 1e-8, rtol = 1e-8)
    sum(sol[, -1]^2)
  }
  obj <- try(MakeADFun(f, list(lr = rep(log(0.3), ns)), silent = TRUE),
             silent = TRUE)
  if (inherits(obj, "try-error")) return(c(NA, NA))
  v <- suppressWarnings(try(obj$fn(obj$par), silent = TRUE))
  g <- suppressWarnings(try(max(abs(obj$gr(obj$par))), silent = TRUE))
  c(if (inherits(v, "try-error")) NA else v,
    if (inherits(g, "try-error")) NA else g)
}

# exact reference for the value
exact <- function(ns) sum(outer(times, rep(0.3, ns), function(t, r)
  exp(-r * t))^2)

for (method in c("lsoda", "lsode", "adams", "rk4", "ode45", "euler")) {
  cat("\n== method:", method, "==\n")
  for (ns in c(1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32)) {
    x <- run(ns, method)
    cat(sprintf("  ns=%3d  fn=%-14s exact=%-12s max|gr|=%s\n", ns,
                format(x[1], digits = 8), format(exact(ns), digits = 8),
                format(x[2], digits = 6)))
    if (is.na(x[1])) break
  }
}
